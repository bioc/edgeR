catchKallisto <- function(parent.dir=NULL,sample.dirs=NULL,DGEList=TRUE,divide=FALSE,verbose=TRUE)
#	Read transcriptwise counts and bootstrap samples from kallisto output
#	Use bootstrap samples to estimate overdispersion of transcriptwise counts
#	Gordon Smyth and Pedro Baldoni
#	Created 2 April 2018. Last modified 29 Aug 2026.
{
#	Check parent.dir
	if(length(parent.dir) > 1L) stop("parent.dir should be of length 1")
	if(is.null(parent.dir)) parent.dir <- "."

#	Check sample.dirs
	if(is.null(sample.dirs)) {
		sample.dirs <- dir(parent.dir)
		IsKallisto <- file.exists(file.path(parent.dir,sample.dirs,"abundance.h5"))
		sample.dirs <- sample.dirs[IsKallisto]
	}

#	Full paths
	paths <- file.path(parent.dir,sample.dirs)
	NSamples <- length(paths)

#	Use rhdf5 package for reading
	suppressPackageStartupMessages(OK <- requireNamespace("rhdf5",quietly=TRUE))
	if(!OK) stop("rhdf5 package required but is not installed (or can't be loaded)")
	
#	Initialize vector of inferential sample types
	ResampleType <- rep_len("bootstrap",NSamples)

#	Accumulate counts and CV^2 of bootstrap counts for each sample
	for (j in 1L:NSamples) {
		if(verbose) cat("Reading ",paths[j],", ",sep="")

#		Open H5 file
		h5File <- file.path(paths[j],"abundance.h5")
		if(!file.exists(h5File)) stop("abundance.h5 file not found at specified path")
		h5 <- rhdf5::H5Fopen(h5File)

#		Auxiliary information
		aux <- h5$aux
		NTx <- length(aux$lengths)
		NBoot <- as.integer(aux$num_bootstrap)
		if(verbose) cat(NTx,"transcripts,",NBoot,"bootstraps\n")

#		Store counts
		if(j == 1L) {
			Counts <- Length <- matrix(0,NTx,NSamples)
			DF <- rep_len(0L,NTx)
			OverDisp <- rep_len(0,NTx)
		}
		Counts[,j] <- h5$est_counts
		Length[,j] <- aux$eff_lengths

#		Bootstraps
		if(NBoot > 0L) Boot <- do.call(cbind,h5$bootstrap)

#		Close H5 file
		rhdf5::H5Fclose(h5)

#		Summarize bootstrap samples
		if(NBoot > 0L) {
			M <- rowMeans(Boot)
			i <- (M > 0)
			OverDisp[i] <- OverDisp[i] + rowSums((Boot[i,]-M[i])^2) / M[i]
			DF[i] <- DF[i]+NBoot-1L
		}
	}
	
# Compute length statistics
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
	Ann <- data.frame(
	    Length=as.integer(aux$lengths),
	    AveTxLength=AveTxLength,
	    Max2MinLength=RangeTxLength,
	    Overdispersion=OverDisp,
	    row.names=aux$ids,
	    stringsAsFactors=FALSE)
	dimnames(Counts) <- list(aux$ids,basename(paths))

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
