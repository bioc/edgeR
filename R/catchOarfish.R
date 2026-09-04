catchOarfish <- function(parent.dir=NULL,prefixes=NULL,DGEList=TRUE,divide=FALSE,verbose=TRUE)
#	Read transcriptwise counts and bootstrap samples from Oarfish output
#	Use bootstrap samples to estimate overdispersion of transcriptwise counts
#	Will unpack Genecode Tx annotation if found in row.names.
#	Gordon Smyth and Pedro Baldoni
#	Created 4 Jul 2025. Last modified 29 Aug 2026.
{
#	Check prefixes
	if(is.null(prefixes)) {
		if(is.null(parent.dir)) {
			QuantFiles <- dir(pattern="*\\.quant$")
		} else {
			QuantFiles <- dir(path=parent.dir,pattern="*\\.quant$")
		}
		n <- nchar(QuantFiles)
		prefixes <- substring(QuantFiles,1,n-6L) 
	} else {
		prefixes <- as.character(prefixes)
	}
	NSamples <- length(prefixes)
	if(NSamples < 1L) stop("No oarfish output files", call.=FALSE)
	if(!is.null(parent.dir)) prefixes <- file.path(parent.dir,prefixes)

#	Use jsonlite and arrow packages for reading
	OK <- requireNamespace("jsonlite",quietly=TRUE)
	if(!OK) stop("jsonlite package required but is not installed (or can't be loaded)")
	OK <- requireNamespace("readr",quietly=TRUE)
	if(!OK) stop("readr package required but is not installed (or can't be loaded)")
	OK <- requireNamespace("nanoparquet",quietly=TRUE)
	if(!OK) stop("nanoparquet package required but is not installed (or can't be loaded)")

#	Initialize vector of inferential sample types
	ResampleType <- rep_len("bootstrap",NSamples)

#	Accumulate counts and CV^2 of bootstrap counts for each sample
	for (j in 1L:NSamples) {
		if(verbose) cat("Reading ",prefixes[j],", ", sep="")

#		File locations
		MetaFile <- paste0(prefixes[j],".meta_info.json")
		QuantFile <- paste0(prefixes[j],".quant")
		BootFile <- paste0(prefixes[j],".infreps.pq")
		if(!file.exists(QuantFile)) stop("quant file not found at specified path")

#		Meta information
		Meta <- jsonlite::fromJSON(MetaFile)
		NBoot <- Meta$num_bootstraps
		if(is.null(NBoot)) stop("Can't find number of bootstraps")
		if(verbose) cat(NBoot,"bootstraps\n")

#		Read counts
		if(j == 1L) {
			Quant <- suppressWarnings(readr::read_tsv(QuantFile,col_types="cdd",progress=FALSE))
			NTx <- nrow(Quant)
			Counts <- matrix(0,NTx,NSamples)
			DF <- rep_len(0L,NTx)
			OverDisp <- rep_len(0,NTx)
			Counts[,1L] <- Quant$num_reads
			Ann <- data.frame(Length=Quant$len)
			row.names(Ann) <- Quant$tname
		} else {
			Quant <- suppressWarnings(readr::read_tsv(QuantFile,col_types="__d",progress=FALSE))
			Counts[,j] <- Quant$num_reads
		}

#		Bootstrap samples
		if(NBoot > 0L) {
			Boot <- as.matrix(nanoparquet::read_parquet(BootFile))
			M <- rowMeans(Boot)
			i <- (M > 0)
			OverDisp[i] <- OverDisp[i] + rowSums((Boot[i,]-M[i])^2) / M[i]
			DF[i] <- DF[i]+NBoot-1L
		}
	}

#	Estimate overdispersion for each transcript
	i <- (DF > 0L)
	if(sum(i) > 0L) {
		OverDisp[i] <- OverDisp[i] / DF[i]
#		Apply a limited amount of moderation
		DFMedian <- median(DF[i])
		DFPrior <- 3
		OverDispPrior <- median(OverDisp[i]) / qf(0.5,df1=DFMedian,df2=DFPrior)
		if(OverDispPrior < 1) OverDispPrior <- 1
		OverDisp[i] <- (DFPrior * OverDispPrior + DF[i]*OverDisp[i]) / (DFPrior + DF[i])
		OverDisp <- pmax(OverDisp,1)
		OverDisp[!i] <- OverDispPrior
	} else {
		OverDisp[] <- NA_real_
		OverDispPrior <- NA_real_
	}

#	Prepare output
	dimnames(Counts) <- list(row.names(Ann),basename(prefixes))
	Ann$Overdispersion <- OverDisp

#	Detect and unpack Gencode tx names
	x <- row.names(Ann)[1]
	gencode <- (nchar(x)-nchar(gsub("|","",x,fixed=TRUE)) >= 8L)
	if(gencode) {
		A <- splitGencodeTxNames(row.names(Ann))
		Ann <- data.frame(Ann,A[,-1])
		row.names(Ann) <- row.names(Counts) <- A[,1]
	}

	#	Divided counts
	if(divide) Counts <- Counts / Ann$Overdispersion
	
	if(DGEList) {
	  y  <- DGEList(count=Counts,genes=Ann)
	  y$overdispersion.prior <- OverDispPrior
	  y$resample.type <- ResampleType
	  y$divided.counts <- divide
	} else {
	  y <- list(counts=Counts,annotation=Ann,overdispersion.prior=OverDispPrior,resample.type=ResampleType,divided.counts=divide)
	}
	
	y

}
