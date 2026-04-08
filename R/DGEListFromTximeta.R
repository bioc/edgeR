DGEListFromTximeta <- function(txm, samples = NULL, group = NULL, remove.zeros = FALSE, divide = FALSE)
# Create DGEList from tximeta() output.
# Pedro Baldoni and Gordon Smyth
# Created 5 Mar 2026. Last modified 8 Apr 2026.
{
	if(!requireNamespace("SummarizedExperiment", quietly = TRUE))
		stop("SummarizedExperiment package required but is not installed (or can't be loaded)")

	if(!is(txm, "SummarizedExperiment"))
		stop("txm object is not of the SummarizedExperiment class")
	
#	Count matrix
	if(!("counts" %in% SummarizedExperiment::assayNames(txm)))
		stop("txm object doesn't contain counts assay")
	counts <- SummarizedExperiment::assay(txm,"counts")
	NSamples <- ncol(counts)
	NGenes <- nrow(counts)

#	Optional filtering of all zero rows
	if(remove.zeros) {
		AllZero <- which(rowMeans(counts) < 1e-6)
		if(length(AllZero)) {
			txm <- txm[-AllZero,]
			counts <- counts[-AllZero,,drop=FALSE]
			NGenes <- nrow(counts)
		}
	}

#	Gene annotation
	if(is(SummarizedExperiment::rowRanges(txm), "GRanges")) {
		genes <- as.data.frame(SummarizedExperiment::rowRanges(txm))
	} else {
		genes <- as.data.frame(SummarizedExperiment::rowData(txm))
	}
	genes$gene_id <- as.character(genes$gene_id)
	
#	Add tx length annotation
	if(!("length" %in% SummarizedExperiment::assayNames(txm)))
		stop("txm object doesn't contain length assay")
	LTxL <- log(SummarizedExperiment::assay(txm, "length"))
	AveTxLength <- exp(rowMeans(LTxL))
	MinLLen <- apply(LTxL, 1, min)
	MaxLLen <- apply(LTxL, 1, max)
	RangeTxLength <- exp(MaxLLen - MinLLen)
	genes$AveLength <- AveTxLength
	genes$Max2MinLength <- RangeTxLength
	
#	RTA overdispersion
	infReps <- grep("infRep",SummarizedExperiment::assayNames(txm))
	if(length(infReps)) {
#		Accumulate means and variances over infReps, using Welford's algorithm
		N <- 1
		M <- SummarizedExperiment::assay(txm,infReps[1],withDimnames=FALSE)
		V <- array(0,dim(M))
		for (j in infReps[-1]) {
			N <- N+1
			Delta <- SummarizedExperiment::assay(txm,j,withDimnames=FALSE)-M
			M <- M + Delta/N
			V <- V + Delta*(SummarizedExperiment::assay(txm,j,withDimnames=FALSE)-M)
		}
		OverDisp <- rowSums(V/M,na.rm=TRUE)
		DF <- (N-1)*rowSums(M>0)
#		Estimate overdispersions
		i <- (DF > 0L)
		if(sum(i) > 0L) {
			OverDisp[i] <- OverDisp[i] / DF[i]
#			Apply a limited amount of moderation
			DFMedian <- median(DF[i])
			DFPrior <- 3
			OverDispPrior <- median(OverDisp[i]) / qf(0.5, df1 = DFMedian, df2 = DFPrior)
			if(OverDispPrior < 1) OverDispPrior <- 1
			OverDisp[i] <- (DFPrior * OverDispPrior + DF[i] * OverDisp[i]) / (DFPrior + DF[i])
			OverDisp <- pmax(OverDisp, 1)
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
			counts <- counts / genes$Overdispersion
		} else {
			divide <- FALSE
		}
	}

#	Construct DGEList
	y <- DGEList(counts=counts,samples=samples,group=group,genes=genes,remove.zeros=FALSE)
	y$divided.counts <- divide

#	Offset matrix
	if(identical(S4Vectors::metadata(txm)$countsFromAbundance, "no")) {
		y$tximport.counts <- "raw"
		PriorOffset <- LTxL - rowMeans(LTxL)
		y$offset.prior <- PriorOffset
	} else {
		y$tximport.counts <- S4Vectors::metadata(txm)$countsFromAbundance
	}

	y
}
