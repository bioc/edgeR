catchSalmon <- function(parent.dir=NULL,sample.dirs=NULL,DGEList=TRUE,divide=FALSE,offset.prior=TRUE,verbose=TRUE)
#	Read transcriptwise counts and bootstrap samples from Salmon output.
#	Use Gibbs or bootstrap samples to estimate overdispersion of transcriptwise counts.
#	Will unpack Genecode Tx annotation if found in row.names.
#	Gordon Smyth and Pedro Baldoni
#	Created 1 April 2018. Last modified 28 Aug 2026.
{
#	Check parent.dir
	if(length(parent.dir) > 1L) stop("parent.dir should be of length 1")
	if(is.null(parent.dir)) parent.dir <- "."

#	Check sample.dirs
	if(is.null(sample.dirs)) {
		sample.dirs <- dir(parent.dir)
		IsSalmon <- file.exists(file.path(parent.dir,sample.dirs,"aux_info"))
		sample.dirs <- sample.dirs[IsSalmon]
	}

#	Full paths
	paths <- file.path(parent.dir,sample.dirs)
	NSamples <- length(paths)

#	Use jsonlite and readr packages for reading
	OK <- requireNamespace("jsonlite",quietly=TRUE)
	if(!OK) stop("jsonlite package required but is not installed (or can't be loaded)")
	OK <- requireNamespace("readr",quietly=TRUE)
	if(!OK) stop("readr package required but is not installed (or can't be loaded)")

#	Initialize vector of inferential sample types
	ResampleType <- rep_len("",NSamples)

#	Accumulate counts and CV^2 of bootstrap counts for each sample
	for (j in 1L:NSamples) {
		if(verbose) cat("Reading ",paths[j],", ",sep="")

#		File locations
		MetaFile <- file.path(paths[j],"aux_info","meta_info.json")
		QuantFile <- file.path(paths[j],"quant.sf")
		BootFile <- file.path(paths[j],"aux_info","bootstrap","bootstraps.gz")
		if(!file.exists(QuantFile)) {
			QuantFile <- dir(paths[j],pattern="^quant.sf",full.names=TRUE)
			if(length(QuantFile)) QuantFile <- QuantFile[1]
			if(!file.exists(QuantFile)) {
				stop("quant.sf file not found at specified path")
			}
		}

#		Meta information
		Meta <- jsonlite::fromJSON(MetaFile)
		NTx <- Meta$num_targets
		if(is.null(NTx)) NTx <- Meta$num_valid_targets
		if(is.null(NTx)) stop("Can't find number of targets")
		NBoot <- Meta$num_bootstraps
		if(is.null(NBoot)) stop("Can't find number of bootstraps")
		Type <- Meta$samp_type
		if(is.null(ResampleType)) Type <- "bootstrap" else ResampleType[j] <- Type
		if(verbose) cat(NTx,"transcripts,",NBoot,Type,"samples\n")

#		Read counts and lengths
		if(j == 1L) {
			Counts <- Length <- matrix(0,NTx,NSamples)
			DF <- rep_len(0L,NTx)
			OverDisp <- rep_len(0,NTx)
			Quant1 <- suppressWarnings(readr::read_tsv(QuantFile,col_types="cdd_d",progress=FALSE))
			Counts[,1L] <- Quant1$NumReads	
			Length[,1L] <- Quant1$EffectiveLength
		} else {
			Quant <- suppressWarnings(readr::read_tsv(QuantFile,col_types="__d_d",progress=FALSE))
			Counts[,j] <- Quant$NumReads
			Length[,j] <- Quant$EffectiveLength
		}

#		Bootstrap samples
		if(NBoot > 0L) {
			BootFileCon <- gzcon(file(BootFile,open="rb"))
			Boot <- readBin(BootFileCon,what="double",n=NTx*NBoot)
			close(BootFileCon)
			dim(Boot) <- c(NTx,NBoot)
			M <- rowMeans(Boot)
			i <- (M > 0)
			OverDisp[i] <- OverDisp[i] + rowSums((Boot[i,]-M[i])^2) / M[i]
			DF[i] <- DF[i]+NBoot-1L
		}
	}
	
#	Compute length statistics
	LTxL <- log(Length)
	AveTxLength <- exp(rowMeans(LTxL))
	MinLLen <- apply(LTxL, 1, min)
	MaxLLen <- apply(LTxL, 1, max)
	RangeTxLength <- exp(MaxLLen - MinLLen)

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
	Quant1 <- as.data.frame(Quant1,stringsAsFactors=FALSE)
	dimnames(Counts) <- list(Quant1$Name,basename(paths))
	row.names(Quant1) <- Quant1$Name
	Quant1$Name <- Quant1$EffectiveLength <- Quant1$NumReads <- NULL
	Quant1$AveLength <- AveTxLength
	Quant1$Max2MinLength <- RangeTxLength
	Quant1$Overdispersion <- OverDisp

#	Detect and unpack Gencode tx names
	x <- row.names(Quant1)[1]
	gencode <- (nchar(x)-nchar(gsub("|","",x,fixed=TRUE)) >= 8L)
	if(gencode) {
		A <- splitGencodeTxNames(row.names(Quant1))
		Quant1 <- data.frame(Quant1,A[,-1])
		row.names(Quant1) <- row.names(Counts) <- A[,1]
	}

#	Divided counts
	if(divide) Counts <- Counts / Quant1$Overdispersion

	if(DGEList) {
		y <- DGEList(count=Counts,genes=Quant1)
		y$overdispersion.prior <- OverDispPrior
		y$resample.type <- ResampleType
		y$divided.counts <- divide
		if(offset.prior) {
			y$offset.prior <- LTxL - rowMeans(LTxL)
			dimnames(y$offset.prior) <- dimnames(Counts)
		}
	} else {
		y <- list(counts=Counts,
			length=Length,
			annotation=Quant1,
			overdispersion.prior=OverDispPrior,
			resample.type=ResampleType,
			divided.counts=divide)
	}
	
	y
}
