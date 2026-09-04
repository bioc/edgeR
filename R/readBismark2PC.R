readBismark2PC <- function(files,sample.names=NULL,readr=TRUE,verbose=TRUE)
#	Read Bismark coverage files and create a PCList
#
#	It is assumed that genomic loci can be represented as integers, so
#	the largest locus position must be less than about 2*10^9.
#	The number of chromosomes times the largest locus position must be
#	less than 10^16.
#
#	Lizhong Chen, Gordon Smyth
#	Created 23 Apr 2025. Last modified 23 Apr 2025.
{
	files <- as.character(files)
	nsamples <- length(files)
	if(is.null(sample.names)) sample.names <- removeExt(removeExt(removeExt(files)))
	if(readr) {
		OK <- requireNamespace("readr",quietly=TRUE)
		if(!OK) stop("readr package required but is not installed (or can't be loaded)")
	}

#	Read all files and store
	CountList1 <- list()
	CountList2 <- list()
	ChrRleList <- list()
	LocusList <- list()
	ChrNames <- c()
	MaxLocus <- 1L
	for(i in seq_len(nsamples)) {
		if(verbose) cat("Reading",files[i],"\n")
		if(readr)
			x <- as.data.frame(suppressWarnings(readr::read_tsv(files[i],col_names=FALSE,col_types="ci__ii",progress=FALSE)))
		else
			x <- read.delim(files[i],header=FALSE,colClasses=c("character","integer","NULL","NULL","integer","integer"))
		ChrRleList[[i]] <- rle(x[,1])
		LocusList[[i]]  <- x[,2]
		CountList1[[i]] <- x[,3]
		CountList2[[i]] <- x[,4]
		ChrNames <- unique(c(ChrNames,ChrRleList[[i]]$values))
	}

	if(verbose) cat("Hashing ...\n")

#	Convert rle values to integer
	for(i in seq_len(nsamples)) ChrRleList[[i]]$values <- match(ChrRleList[[i]]$values,ChrNames)

#	Hash the genomic positions
	HashBase <- length(ChrNames)+1L
	HashList <- list()
	HashUnique <- c()
	for (i in seq_len(nsamples)) {
		HashList[[i]] <- inverse.rle(ChrRleList[[i]]) / HashBase + LocusList[[i]]
		HashUnique <- unique(c(HashUnique,HashList[[i]]))
	}

	if(verbose) cat("Collating counts ...\n")

#	Merged count matrix
	counts1 <- counts2 <- matrix(0L,length(HashUnique),nsamples)
	for(i in seq_len(nsamples)) {
		m <- match(HashList[[i]], HashUnique)
		counts1[m,i] <- CountList1[[i]]
		counts2[m,i] <- CountList2[[i]]
	}

#	Unhash
	Locus <- as.integer(HashUnique)
	Chr <- as.integer( (HashUnique-Locus) * HashBase + 0.5 )
	attr(Chr,"levels") <- ChrNames
	class(Chr) <- "factor"

#	Attach dimension names and form DGEList
	colnames(counts1) <- colnames(counts2) <- sample.names
	y <- PCList(counts1, counts2, genes=data.frame(Chr,Locus))
	row.names(y) <- paste(Chr,Locus,sep="-")
	y
}

