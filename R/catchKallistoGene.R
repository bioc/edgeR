catchKallistoGene <- function(parent.dir=NULL,sample.dirs=NULL,tx2gene=NULL,remove.version.numbers=TRUE,DGEList=TRUE,divide=FALSE,impute.eff.len=TRUE,offset.prior=TRUE,gene.length="moderate",verbose=TRUE)
#	Read transcriptwise counts and bootstrap samples from kallisto output
#	and summarize at gene level using either imbedded Gencode annotation
#	or an externally provided data.frame mapping tx to gene IDs.
#	Use bootstrap samples to estimate overdispersion of genewise counts.
#	Gordon Smyth
#	Created 11 Sep 2026. Last modified 12 Sep 2026.
{
#	Check specified directories
	if(length(parent.dir) > 1L) stop("parent.dir should be of length 1")
	if(is.null(sample.dirs)) {
		if(is.null(parent.dir)) parent.dir <- "."
		sample.dirs <- dir(parent.dir)
		IsKallisto <- file.exists(file.path(parent.dir,sample.dirs,"abundance.h5"))
		sample.dirs <- sample.dirs[IsKallisto]
	}

#	Full paths
	if(is.null(parent.dir)) {
		paths <- sample.dirs
	} else {
		paths <- file.path(parent.dir,sample.dirs)
	}
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
		NBoot <- as.integer(aux$num_bootstrap)
		if(verbose) cat(length(aux$ids),"transcripts,",NBoot,"bootstraps\n")

#		Initalize dataset-wide information from first sample
		if(j == 1L) {
			TxID <- aux$ids
			NTx <- length(TxID)
			TxLen <- as.vector(aux$lengths)

#			Get gene IDs from tx2gene or from unpacking Gencode annotation
			TxID1 <- TxID[1]
			if(is.null(tx2gene)) {
				IsGencode <- (nchar(TxID1)-nchar(gsub("|","",TxID1,fixed=TRUE)) >= 8L)
				if(IsGencode)
					if(verbose) message("Summarizing to genewise counts using Gencode's imbedded annotation")
				else
					stop("tx2gene not provided and row names do not appear to be from Gencode")
				GenecodeAnn <- splitGencodeTxNames(TxID, remove.version.numbers=remove.version.numbers)
				EnsG <- GenecodeAnn[,"EnsG"]
				d <- duplicated(EnsG)
				GeneAnn <- data.frame(GenecodeAnn[!d,"GeneName",drop=FALSE])
			} else {
				tx2gene <- as.data.frame(tx2gene)
				if(ncol(tx2gene) < 2L) stop("tx2gene doesn't have two columns")
				if(remove.version.numbers) {
					tx2gene[,1] <- strsplit2(tx2gene[,1],split="\\.")[,1]
					TxID <- strsplit2(TxID,split="\\.")[,1]
				}
				m <- match(TxID,tx2gene[,1])
				if(anyNA(m)) stop("Tx names not found in first column of tx2gene")
				EnsG <- tx2gene[m,2]
				if(anyNA(EnsG)) stop("Missing gene IDs")
				GeneAnn <- NULL	
				if(verbose) message("Summarizing to genewise counts using tx2gene")
			}

#			Genewise statistics
			NTxPerGene <- rowsum(rep_len(1L,NTx),EnsG,reorder=FALSE)
			EnsGu <- row.names(NTxPerGene)
			NGene <- length(NTxPerGene)
			NTxPerGene <- drop(NTxPerGene)

#			Initialize matrices
			Counts <- matrix(0,NGene,NSamples)
			TPM <- EffTxLen <- matrix(0,NTx,NSamples)
			DF <- rep_len(0L,NGene)
			OverDisp <- rep_len(0,NGene)
		}

#		Store quantifications
		Counts[,j] <- drop(rowsum(h5$est_counts,EnsG,reorder=FALSE))
		EffTxLen[,j] <- aux$eff_lengths
		x <- h5$est_counts / aux$eff_lengths
		TPM[,j] <- 1e6 * x / sum(x)

#		Bootstrap samples
		if(NBoot > 0L) {
			Boot <- do.call(cbind,h5$bootstrap)
			Boot <- rowsum(Boot,EnsG,reorder=FALSE)
			M <- rowMeans(Boot)
			i <- (M > 0)
			OverDisp[i] <- OverDisp[i] + rowSums((Boot[i,]-M[i])^2) / M[i]
			DF[i] <- DF[i]+NBoot-1L
		}

#		Close H5 file
		rhdf5::H5Fclose(h5)
	}

#	Maximum tx length per gene
	o <- order(TxLen,decreasing=TRUE)
	m <- match(EnsGu,EnsG[o])
	MaxTxLen <- TxLen[o][m]

#	Impute effective lengths
	if(impute.eff.len) EffTxLen <- .imputeEffectiveLengths(TxLen,EffTxLen)

#	Average gene length, with weak moderation towards genewise average and towards unweighted average
	gene.length <- match.arg(gene.length,c("moderate","tximport","simple"))
	if(identical(gene.length,"moderate")) {
		eps <- 1e-6
		m <- rowMeans(TPM)
		TPM2 <- TPM + eps + m/100
		EffGeneLen <- rowsum(TPM2*EffTxLen,EnsG,reorder=FALSE)/rowsum(TPM2,EnsG,reorder=FALSE)
	}
	if(identical(gene.length,"tximport")) {
		GeneIsAllZero <- which(rowSums(Counts) == 0)
		TxGeneIsAllZero <- which(EnsG %in% EnsG[!d][GeneIsAllZero])
		m <- rowMeans(EffTxLen[TxGeneIsAllZero,,drop=FALSE])
		EffGeneLen <- Counts
		EffGeneLen[GeneIsAllZero,] <- rowsum(m,EnsG[TxGeneIsAllZero],reorder=FALSE) / NTxPerGene[GeneIsAllZero]
		EffGeneLen[-GeneIsAllZero,] <- rowsum(TPM[-TxGeneIsAllZero,]*EffTxLen[-TxGeneIsAllZero,],EnsG[-TxGeneIsAllZero],reorder=FALSE) / rowsum(TPM[-TxGeneIsAllZero,],EnsG[-TxGeneIsAllZero],reorder=FALSE) 
		if(anyNA(EffGeneLen)) {
			m <- exp(rowMeans(log(EffGeneLen),na.rm=TRUE))
			i <- which(is.na(EffGeneLen))
			EffGeneLen[i] <- matrix(m,NGene,NSamples)[i]
		}
	}
	if(identical(gene.length,"simple")) {
		EffGeneLen <- rowsum(EffTxLen,EnsG,reorder=FALSE)/ NTxPerGene
	}

#	Compute length statistics
	LGL <- log(EffGeneLen)
	m <- rowMeans(LGL)
	AveLength <- exp(m)
	MinLLen <- apply(LGL, 1, min)
	MaxLLen <- apply(LGL, 1, max)
	RangeLength <- exp(MaxLLen - MinLLen)

#	Estimate overdispersion for each transcript or gene
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
	dimnames(Counts) <- dimnames(EffGeneLen) <- list(EnsGu,basename(paths))
	NTxPerGene <- rowsum(rep_len(1L,NTx),EnsG,reorder=FALSE)
	if(is.null(GeneAnn))
		Genes <- data.frame(NTx=NTxPerGene,MaxTxLen=MaxTxLen,AveEffLen=AveLength,Max2MinEffLen=RangeLength,Overdispersion=OverDisp)
	else
		Genes <- data.frame(GeneAnn,NTx=NTxPerGene,MaxTxLen=MaxTxLen,AveEffLen=AveLength,Max2MinEffLen=RangeLength,Overdispersion=OverDisp)
	row.names(Genes) <- EnsGu

#	Divided counts
	if(divide) {
		if(NBoot > 0) {
			Counts <- Counts / OverDisp
		} else {
			message("No bootscript or Gibbs samples, so counts not divided")
		}
	}

	if(DGEList) {
		y <- DGEList(count=Counts,genes=Genes)
		y$overdispersion.prior <- OverDispPrior
		y$resample.type <- ResampleType
		y$divided.counts <- divide
		if(offset.prior) {
			y$offset.prior <- LGL - m
			dimnames(y$offset.prior) <- dimnames(Counts)
		}
		y$other$effective.length <- EffGeneLen
	} else {
		y <- list(counts=Counts,
			effective.length=EffGeneLen,
			annotation=Genes,
			overdispersion.prior=OverDispPrior,
			resample.type=ResampleType,
			divided.counts=divide)
	}
	
	y
}
