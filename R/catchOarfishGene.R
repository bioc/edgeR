catchOarfishGene <- function(parent.dir=NULL,prefixes=NULL,tx2gene=NULL,remove.version=TRUE,DGEList=TRUE,divide=TRUE,gene.length="moderate",verbose=TRUE)
#	Read transcriptwise counts and bootstrap samples from Oarfish output
#	and summarize at gene level using either imbedded Gencode annotation
#	or an externally provided data.frame mapping tx to gene IDs.
#	Use bootstrap samples to estimate overdispersion of genewise counts.
#	Gordon Smyth and Pedro Baldoni
#	catchOarfish() created 4 Jul 2025.
#	catchOarfishGene() created 27 Oct 2026. Last modified 31 Oct 2026.
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

#	Use jsonlite and nanoparquet packages for reading
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

#		Read counts and lengths
		if(j == 1L) {
			Quant <- suppressWarnings(readr::read_tsv(QuantFile,col_types="cdd",progress=FALSE))
			NTx <- nrow(Quant)
#			Detect and unpack Gencode tx names
			Name1 <- Quant$tname[1]
			IsGencode <- (nchar(Name1)-nchar(gsub("|","",Name1,fixed=TRUE)) >= 8L)
			if(is.null(tx2gene)) {
				if(!IsGencode) stop("tx2gene not provided and row names do not appear to be from Gencode")
				if(verbose) message("Summarizing to genewise counts using Gencode's imbedded annotation")
				TxAnn <- splitGencodeTxNames(Quant$tname, remove.version=remove.version)
				EnsG <- TxAnn[,"EnsG"]
				GeneAnn <- collapseGencodeTx2Gene(TxAnn)
			} else {
				tx2gene <- as.data.frame(tx2gene)
				if(ncol(tx2gene) < 2L) stop("tx2gene doesn't have two columns")
				if(IsGencode) {
					tx2geneIsGencode <- (nchar(tx2gene[1,1])-nchar(gsub("|","",tx2gene[1,1],fixed=TRUE)) >= 8L)
					if(IsGencode && !tx2geneIsGencode) {
						if(verbose) message("Salmon tx names are Gencode-style but tx2gene tx names are not, so unpacking Salmon tx names")
						Quant$tname <- sub("[|].*", "", Quant$tname)
					}
				}
				if(remove.version) {
					tx2gene[,1] <- sub("\\.[0-9]+", "", tx2gene[,1])
					tx2gene[,2] <- sub("\\.[0-9]+", "", tx2gene[,2])
					Quant$tname <- sub("\\.[0-9]+", "", Quant$tname)
				}
				m <- match(Quant$tname,tx2gene[,1])
				if(anyNA(m)) stop("Tx names not found in first column of tx2gene")
				EnsG <- tx2gene[m,2]
				if(anyNA(EnsG)) stop("Missing gene IDs")
				GeneAnn <- NULL
				if(verbose) message("Summarizing to genewise counts using tx2gene")
			}
			NTxPerGene <- rowsum(rep_len(1L,NTx),EnsG,reorder=FALSE)
			EnsGu <- row.names(NTxPerGene)
			NGene <- length(NTxPerGene)
			NTxPerGene <- drop(NTxPerGene)
			TxCounts <- matrix(0,NTx,NSamples)
			GeneCounts <- matrix(0,NGene,NSamples)
			DF <- rep_len(0L,NGene)
			OverDisp <- rep_len(0,NGene)
			TxCounts[,1L] <- Quant$num_reads
			GeneCounts[,1L] <- drop(rowsum(Quant$num_reads,EnsG,reorder=FALSE))
			TxLen <- Quant$len
		} else {
			Quant <- suppressWarnings(readr::read_tsv(QuantFile,col_types="__d",progress=FALSE))
			TxCounts[,j] <- Quant$num_reads
			GeneCounts[,j] <- drop(rowsum(Quant$num_reads,EnsG,reorder=FALSE))
		}

#		Bootstrap samples
		if(NBoot > 0L) {
			Boot <- as.matrix(nanoparquet::read_parquet(BootFile))
			Boot <- rowsum(Boot,EnsG,reorder=FALSE)
			M <- rowMeans(Boot)
			i <- (M > 0)
			OverDisp[i] <- OverDisp[i] + rowSums((Boot[i,]-M[i])^2) / M[i]
			DF[i] <- DF[i]+NBoot-1L
		}
	}

#	Maximum tx length per gene
	o <- order(TxLen,decreasing=TRUE)
	m <- match(EnsGu,EnsG[o])
	MaxTxLen <- TxLen[o][m]

#	Compute TPM
	TPM <- t(t(TxCounts) / colSums(TxCounts)) * 1e6
	TPMGene <- rowsum(TPM,EnsG,reorder=FALSE)

#	Average gene length, with weak moderation towards genewise average and towards unweighted average
	gene.length <- match.arg(gene.length,c("moderate","tximport","simple"))
	if(identical(gene.length,"moderate")) {
		eps <- 1e-6
		m <- rowMeans(TPM)
		TPM2 <- TPM + eps + m/100
		GeneLen <- rowsum(TPM2*TxLen,EnsG,reorder=FALSE)/rowsum(TPM2,EnsG,reorder=FALSE)
	}
	if(identical(gene.length,"tximport")) {
		GeneIsAllZero <- which(rowSums(GeneCounts) == 0)
		TxGeneIsAllZero <- which(EnsG %in% EnsGu[GeneIsAllZero])
		m <- TxLen[TxGeneIsAllZero]
		GeneLen <- GeneCounts
		GeneLen[GeneIsAllZero,] <- c(rowsum(m,EnsG[TxGeneIsAllZero],reorder=FALSE) / NTxPerGene[GeneIsAllZero])
		GeneLen[-GeneIsAllZero,] <- rowsum(TxCounts[-TxGeneIsAllZero,]*TxLen[-TxGeneIsAllZero],EnsG[-TxGeneIsAllZero],reorder=FALSE) / rowsum(TxCounts[-TxGeneIsAllZero,],EnsG[-TxGeneIsAllZero],reorder=FALSE)
		if(anyNA(GeneLen)) {
			m <- exp(rowMeans(log(GeneLen),na.rm=TRUE))
			i <- which(is.na(GeneLen))
			GeneLen[i] <- matrix(m,NGene,NSamples)[i]
		}
	}
	if(identical(gene.length,"simple")) {
		GeneLen <- rowsum(TxLen,EnsG,reorder=FALSE)/ NTxPerGene
	}

#	Compute length statistics
	LGL <- log(GeneLen)
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
	dimnames(GeneCounts) <- dimnames(GeneLen) <- dimnames(TPMGene) <- list(EnsGu,basename(prefixes))
	NTxPerGene <- rowsum(rep_len(1L,NTx),EnsG,reorder=FALSE)
	if(is.null(GeneAnn))
		Genes <- data.frame(NTx=NTxPerGene,MaxTxLen=MaxTxLen,AveEffLen=AveLength,Overdispersion=OverDisp)
	else
		Genes <- data.frame(GeneAnn,NTx=NTxPerGene,MaxTxLen=MaxTxLen,AveEffLen=AveLength,Overdispersion=OverDisp)
	row.names(Genes) <- EnsGu

#	Divided counts
	if(divide) {
		if(NBoot > 0) {
			GeneCounts <- GeneCounts / OverDisp
		} else {
			message("No bootstrap samples, so counts not divided")
		}
	}

	if(DGEList) {
		y <- DGEList(count=GeneCounts,genes=Genes)
		y$overdispersion.prior <- OverDispPrior
		y$resample.type <- ResampleType
		y$divided.counts <- divide
		y$other$effective.length <- GeneLen
		y$other$tpm <- TPMGene
	} else {
		y <- list(counts=GeneCounts,
			effective.length=GeneLen,
			tpm=TPMGene,
			annotation=Genes,
			overdispersion.prior=OverDispPrior,
			resample.type=ResampleType,
			divided.counts=divide)
	}

	y
}
