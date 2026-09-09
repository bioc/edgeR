catchSalmonGene <- function(parent.dir=NULL,sample.dirs=NULL,tx2gene=NULL,remove.version.numbers=TRUE,DGEList=TRUE,divide=FALSE,impute.eff.len=TRUE,offset.prior=TRUE,gene.length="moderate",verbose=TRUE)
#	Read transcriptwise counts and bootstrap samples from Salmon output
#	and summarize at gene level using either imbedded Gencode annotation
#	or an externally provided data.frame mapping tx to gene IDs.
#	Use Gibbs or bootstrap samples to estimate overdispersion of genewise counts.
#	Gordon Smyth and Pedro Baldoni
#	catchSalmon() created 1 Apr 2018. 
#	catchSalmonWithGencode() created 14 July 2026.
#	catchSalmonWithGene() created 6 Sep 2026. Last modified 7 Sep 2026.
{
#	Check specified directories
	if(length(parent.dir) > 1L) stop("parent.dir should be of length 1")
	if(is.null(sample.dirs)) {
		if(is.null(parent.dir)) parent.dir <- "."
		sample.dirs <- dir(parent.dir)
		IsSalmon <- file.exists(file.path(parent.dir,sample.dirs,"aux_info"))
		sample.dirs <- sample.dirs[IsSalmon]
	}
	if(is.null(parent.dir)) {
		paths <- sample.dirs
	} else {
		paths <- file.path(parent.dir,sample.dirs)
	}

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
			Quant1 <- suppressWarnings(readr::read_tsv(QuantFile,col_types="cdddd",progress=FALSE))
#			Detect and unpack Gencode tx names
			Name1 <- Quant1$Name[1]
			if(is.null(tx2gene)) {
				IsGencode <- (nchar(Name1)-nchar(gsub("|","",Name1,fixed=TRUE)) >= 8L)
				if(IsGencode)
					if(verbose) message("Summarizing to genewise counts using Gencode's imbedded annotation")
				else
					stop("tx2gene not provided and row names do not appear to be from Gencode")
				GenecodeAnn <- splitGencodeTxNames(Quant1$Name, remove.version.numbers=remove.version.numbers)
				EnsG <- GenecodeAnn[,"EnsG"]
				d <- duplicated(EnsG)
				GeneAnn <- data.frame(GenecodeAnn[!d,"GeneName",drop=FALSE])
			} else {
				tx2gene <- as.data.frame(tx2gene)
				if(ncol(tx2gene) < 2L) stop("tx2gene doesn't have two columns")
				if(remove.version.numbers) {
					tx2gene[,1] <- strsplit2(tx2gene[,1],split="\\.")[,1]
					Quant1$Name <- strsplit2(Quant1$Name,split="\\.")[,1]
				}
				m <- match(Quant1$Name,tx2gene[,1])
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
			Counts <- matrix(0,NGene,NSamples)
			TPM <- EffLen <- matrix(0,NTx,NSamples)
			DF <- rep_len(0L,NGene)
			OverDisp <- rep_len(0,NGene)
			Counts[,1L] <- drop(rowsum(Quant1$NumReads,EnsG,reorder=FALSE))
			TPM[,1L] <- Quant1$TPM
			EffLen[,1L] <- Quant1$EffectiveLength
		} else {
			Quant <- suppressWarnings(readr::read_tsv(QuantFile,col_types="__ddd",progress=FALSE))
			Counts[,j] <- drop(rowsum(Quant$NumReads,EnsG,reorder=FALSE))
			TPM[,j] <- Quant$TPM
			EffLen[,j] <- Quant$EffectiveLength
		}

#		Bootstrap samples
		if(NBoot > 0L) {
			BootFileCon <- gzcon(file(BootFile,open="rb"))
			Boot <- readBin(BootFileCon,what="double",n=NTx*NBoot)
			close(BootFileCon)
			dim(Boot) <- c(NTx,NBoot)
			Boot <- rowsum(Boot,EnsG,reorder=FALSE)
			M <- rowMeans(Boot)
			i <- (M > 0)
			OverDisp[i] <- OverDisp[i] + rowSums((Boot[i,]-M[i])^2) / M[i]
			DF[i] <- DF[i]+NBoot-1L
		}
	}

#	Maximum tx length per gene
	o <- order(Quant1$Length,decreasing=TRUE)
	m <- match(EnsGu,EnsG[o])
	MaxTxLen <- Quant1$Length[o][m]

#	Impute effective lengths
	if(impute.eff.len) EffLen <- .imputeEffectiveLengths(Quant1$Length,EffLen)

#	Average gene length, with weak moderation towards genewise average and towards unweighted average
	gene.length <- match.arg(gene.length,c("moderate","tximport","simple"))
	if(identical(gene.length,"moderate")) {
		eps <- 1e-6
		m <- rowMeans(TPM)
		TPM2 <- TPM + eps + m/100
		Length <- rowsum(TPM2*EffLen,EnsG,reorder=FALSE)/rowsum(TPM2,EnsG,reorder=FALSE)
	}
	if(identical(gene.length,"tximport")) {
		GeneIsAllZero <- which(rowSums(Counts) == 0)
		TxGeneIsAllZero <- which(EnsG %in% EnsG[!d][GeneIsAllZero])
		m <- rowMeans(EffLen[TxGeneIsAllZero,,drop=FALSE])
		Length <- Counts
		Length[GeneIsAllZero,] <- rowsum(m,EnsG[TxGeneIsAllZero],reorder=FALSE) / NTxPerGene[GeneIsAllZero]
		Length[-GeneIsAllZero,] <- rowsum(TPM[-TxGeneIsAllZero,]*EffLen[-TxGeneIsAllZero,],EnsG[-TxGeneIsAllZero],reorder=FALSE) / rowsum(TPM[-TxGeneIsAllZero,],EnsG[-TxGeneIsAllZero],reorder=FALSE) 
		if(anyNA(Length)) {
			m <- exp(rowMeans(log(Length),na.rm=TRUE))
			i <- which(is.na(Length))
			Length[i] <- matrix(m,NGene,NSamples)[i]
		}
	}
	if(identical(gene.length,"simple")) {
		Length <- rowsum(EffLen,EnsG,reorder=FALSE)/ NTxPerGene
	}

#	Compute length statistics
	LGL <- log(Length)
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
	dimnames(Counts) <- dimnames(Length) <- list(EnsGu,basename(paths))
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
		y$other$effective.length <- Length
	} else {
		y <- list(counts=Counts,
			effective.length=Length,
			annotation=Genes,
			overdispersion.prior=OverDispPrior,
			resample.type=ResampleType,
			divided.counts=divide)
	}
	
	y
}
