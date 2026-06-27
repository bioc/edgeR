DGEListFromTximport <- function(txi, samples = NULL, group = NULL, genes = NULL, remove.zeros = FALSE, divide = FALSE)
# Create DGEList from tximport() output.
# Created 2 Feb 2026. Last modified 11 May 2026.
{
#	Check input
	ExpectedCols <- c("counts","length","countsFromAbundance")
	k <- hasName(txi, ExpectedCols)
	if(!all(k))
		stop("Component(s) ",paste(ExpectedCols[!k],collapse=",")," not found in txi")

#	Count matrix
	NSamples <- ncol(txi$counts)
	NGenes <- nrow(txi$counts)

#	Add tx length annotation
	LTxL <- log(txi$length)
	AveLength <- exp(rowMeans(LTxL))
	MinLLen <- apply(LTxL,1,min)
	MaxLLen <- apply(LTxL,1,max)
	Max2MinLength <- exp(MaxLLen - MinLLen)
	if(is.null(genes)) {
		genes <- data.frame(AveLength,Max2MinLength)
	} else {
		genes <- as.data.frame(genes)
		if(!identical(nrow(genes),NGenes)) stop("nrow(genes) is different from nrow(txi$counts)")
		genes$AveLength <- AveLength
		genes$Max2MinLength <- Max2MinLength
	}

#	RTA overdispersion
	if(hasName(txi,"infReps")) {
#		Accumulate genewise overdispersions
		OverDisp <- rep_len(0,NGenes)
		DF <- rep_len(0,NGenes)
		for (j in 1L:NSamples) {
			Boot <- txi$infReps[[j]]
			NBoot <- ncol(Boot)
			M <- rowMeans(Boot)
			i <- (M > 0)
			OverDisp[i] <- OverDisp[i] + rowSums((Boot[i,]-M[i])^2) / M[i]
			DF[i] <- DF[i]+NBoot-1L
		}
#		Estimate overdispersions
		i <- (DF > 0L)
		if(sum(i) > 0L) {
			OverDisp[i] <- OverDisp[i] / DF[i]
			# Apply a limited amount of moderation
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
		genes$Overdispersion <- OverDisp
	}

#	Divided counts
	if(divide) {
		if(hasName(genes,"Overdispersion")) {
			txi$counts <- txi$counts / genes$Overdispersion
		} else {
			divide <- FALSE
		}
	}

#	Construct DGEList
	y <- DGEList(counts=txi$counts,samples=samples,group=group,genes=genes,remove.zeros=FALSE)
	y$divided.counts <- divide

#	Offset matrix
	if(identical(txi$countsFromAbundance,"no")) {
		y$tximport.counts <- "raw"
		y$offset.prior <- LTxL - rowMeans(LTxL)
	} else {
		y$tximport.counts <- txi$countsFromAbundance
	}

	if(remove.zeros) {
		AllZero <- which(rowMeans(y$counts) < 1e-6)
		if(length(AllZero)) y <- y[-AllZero,]
	}

	y
}
