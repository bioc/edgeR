DGEListFromTximport <- function(txi, samples = NULL, group = NULL, genes = NULL, remove.zeros = FALSE, norm.method = "none")
# Create a DGEList object from genewise output from tximport().
# Created 2 Feb 2026. Last modified 3 Feb 2026.
{
#	Count matrix
	NSamples <- ncol(txi$counts)
	NGenes <- nrow(txi$counts)

#	Add tx length annotation
	LTxL <- log(txi$length)
	AveTxLength <- exp(rowMeans(LTxL))
	MinLLen <- apply(LTxL,1,min)
	MaxLLen <- apply(LTxL,1,max)
	RangeTxLength <- exp(MaxLLen - MinLLen)
	if(!is.null(genes)) {
		genes <- as.data.frame(genes)
		if(!identical(nrow(genes),NGenes)) stop("nrow(genes) is different from nrow(txi$counts)")
		genes$TxLengthAve <- AveTxLength
		genes$TxLengthRange <- RangeTxLength
	} else {
		genes <- data.frame(AveTxLength,RangeTxLength)
	}

#	Raw library sizes
	LibSize <- colSums(txi$counts)

#	Normalize library sizes
	NormFactors <- normLibSizes(txi$counts, method=norm.method)

#	RTA overdispersion
	if(!is.null(txi$infReps)) {
#		Accumulate genewise CVs
		OverDisp <- rep_len(1,NGenes)
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

#	Construct DGEList
	y <- DGEList(counts=txi$counts,samples=samples,norm.factors=NormFactors,group=group,genes=genes,remove.zeros=FALSE)

#	Offset matrix
	if(identical(txi$countsFromAbundance,"no")) {
		PriorOffset <- LTxL - rowMeans(LTxL)
		y$offset.prior <- PriorOffset
	}

	if(remove.zeros) {
		AllZero <- which(rowMeans(y$counts) < 1e-6)
		if(length(AllZero)) y <- y[-AllZero,]
	}

	y
}
