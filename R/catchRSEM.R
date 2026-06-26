catchRSEM <- function(files=NULL,path=".",ngibbs=100,DGEList=FALSE,divide=FALSE,verbose=TRUE)
# Read transcriptwise counts and Gibbs posterior means and standard deviations from RSEM output.
# Use Gibbs samples to estimate overdispersion of transcriptwise counts.
# Pedro Baldoni and Gordon Smyth
# Created 24 April 2024. Last modified 26 Jun 2026.
{
#	Check files
	if(is.null(files)) {
		files <- dir(path=path,pattern="*\\.isoforms\\.results$")
	} else {
		files <- as.character(files)
	}
	NSamples <- length(files)
	if(NSamples < 1L) stop("No isoforms.results files", call.=FALSE)
	files <- file.path(path,files)

#	Check ngibbs
	ngibbs <- rep_len(ngibbs,NSamples)
	
#	Initialize vector of inferential sample types
	ResampleType <- rep_len("gibbs",NSamples)

#	Check for readr package  
	OK <- requireNamespace("readr",quietly=TRUE)
	if(!OK) stop("readr package required but is not installed (or can't be loaded)")
	
#	Accumulate counts and CV^2 of bootstrap counts for each sample
	for (j in 1L:NSamples) {
		if(verbose) cat("Reading ",files[j],", ",sep="")
		
#		File location
		QuantFile <- files[j]
		if(!file.exists(QuantFile)) stop("isoforms.results file not found at specified path", call.=FALSE)

#		Read counts
		if(j == 1L) {
			Quant1 <- suppressWarnings(readr::read_tsv(QuantFile,col_types="c_ddd___dd___",progress=FALSE))
			NTx <- nrow(Quant1)
			Counts <- Length <- matrix(0,NTx,NSamples)
			DF <- rep_len(0L,NTx)
			OverDisp <- rep_len(0,NTx)
			if(is.null(Quant1$expected_count)) stop("File doesn't contain expected_count column", call.=FALSE)
			Counts[,1L] <- Quant1$expected_count
			Length[,1L] <- Quant1$effective_length
			M <- Quant1$posterior_mean_count
			S <- Quant1$posterior_standard_deviation_of_count
		} else {
			Quant <- suppressWarnings(readr::read_tsv(QuantFile,col_types="___dd___dd___",progress=FALSE))
			if(is.null(Quant$expected_count)) stop("File doesn't contain expected_count column", call.=FALSE)
			Counts[,j] <- Quant$expected_count
			Length[,j] <- Quant1$effective_length
			M <- Quant$posterior_mean_count
			S <- Quant$posterior_standard_deviation_of_count
		}
		
		NGibbs <- ngibbs[j]
		if(verbose) cat(NTx,"transcripts,",NGibbs,"Gibbs samples\n")
		
#		Bootstrap samples
		NGibbs <- ngibbs[j]
		if(is.null(M)) NGibbs <- 0L
		if(NGibbs > 0L) {
			i <- (M > 0)
			OverDisp[i] <- OverDisp[i] + (NGibbs-1)*(S[i]^2)/M[i]
			DF[i] <- DF[i] + NGibbs-1L
		}
	}
	
# Compute length statistics
	LTxL <- log(Length)
	LTxL[is.infinite(LTxL)] <- log(1e-8)
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
	SampleNames <- removeExt(removeExt(files))
	Quant1 <- as.data.frame(Quant1,stringsAsFactors=FALSE)
	dimnames(Counts) <- list(Quant1$transcript_id,SampleNames)
	row.names(Quant1) <- Quant1$transcript_id
	Quant1$transcript_id <- Quant1$effective_length <- Quant1$expected_count <- NULL
	Quant1$posterior_mean_count <- Quant1$posterior_standard_deviation_of_count<- NULL
	colnames(Quant1) <- c("Length")
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
	  y  <- DGEList(count=Counts,genes=Quant1)
	  y$overdispersion.prior <- OverDispPrior
	  y$resample.type <- ResampleType
	  y$divided.counts <- divide
	} else {
	  y <- list(counts=Counts,annotation=Quant1,overdispersion.prior=OverDispPrior,resample.type=ResampleType,divided.counts=divide)
	}
	
	y
}
