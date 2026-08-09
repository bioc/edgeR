catchSalmonWithGencode <- function(parent.dir=NULL,sample.dirs=NULL,DGEList=TRUE,divide=FALSE,offset.prior=TRUE,verbose=TRUE)
#	Read transcriptwise counts and bootstrap samples from Salmon output.
#	Unpack Gencode annotation and summarize to gene level.
#	Use Gibbs or bootstrap samples to estimate overdispersion of genewise counts.
#	Gordon Smyth and Pedro Baldoni
#	Created 1 April 2018. Last modified 10 Aug 2026.
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
	if(verbose) message("Summarizing to genewise counts using Gencode's imbedded annotation")

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
			IsGencode <- (nchar(Name1)-nchar(gsub("|","",Name1,fixed=TRUE)) >= 8L)
			if(!IsGencode) stop("Row names do not appear to be fromm Gencode")
			GenecodeAnn <- splitGencodeTxNames(Quant1$Name)
			EnsG <- GenecodeAnn[,"EnsG"]
			d <- duplicated(EnsG)
			GeneAnn <- data.frame(GenecodeAnn[!d,c("GeneName","Type")])
			NGene <- nrow(GeneAnn)
			Counts <- Length <- matrix(0,NGene,NSamples)
			DF <- rep_len(0L,NGene)
			OverDisp <- rep_len(0,NGene)
			Counts[,1L] <- drop(rowsum(Quant1$NumReads,EnsG,reorder=FALSE))
			eps <- 1e-6
			Length[,1L] <- rowsum(Quant1$EffectiveLength*(Quant1$TPM+eps),EnsG,reorder=FALSE)/rowsum(Quant1$TPM+eps,EnsG,reorder=FALSE)
		} else {
			Quant <- suppressWarnings(readr::read_tsv(QuantFile,col_types="__ddd",progress=FALSE))
			Counts[,j] <- drop(rowsum(Quant$NumReads,EnsG,reorder=FALSE))
			Length[,j] <- rowsum(Quant$EffectiveLength*(Quant$TPM+eps),EnsG,reorder=FALSE)/rowsum(Quant$TPM+eps,EnsG,reorder=FALSE)
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
	
#	Compute length statistics
	LTxL <- log(Length)
	AveTxLength <- exp(rowMeans(LTxL))
	MinLLen <- apply(LTxL, 1, min)
	MaxLLen <- apply(LTxL, 1, max)
	RangeTxLength <- exp(MaxLLen - MinLLen)

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
	EnsGu <- EnsG[!d]
	dimnames(Counts) <- list(EnsGu,basename(paths))
	NTxPerGene <- rowsum(rep_len(1L,NTx),EnsG,reorder=FALSE)
	Genes <- data.frame(GeneAnn,NTx=NTxPerGene,AveLength=AveTxLength,Max2MinLength=RangeTxLength,Overdispersion=OverDisp)
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
			y$offset.prior <- LTxL - rowMeans(LTxL)
			dimnames(y$offset.prior) <- dimnames(Counts)
		}
	} else {
		y <- list(counts=Counts,
			length=Length,
			annotation=Genes,
			overdispersion.prior=OverDispPrior,
			resample.type=ResampleType,
			divided.counts=divide)
	}
	
	y
}
