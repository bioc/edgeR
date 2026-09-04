rpkm <- function(y, ...)
UseMethod("rpkm")

rpkm.DGEList <- function(y, gene.length=NULL, normalized.lib.sizes=TRUE, log=FALSE, prior.count=2, ...)
#	RPKM for a DGEList.
#	Gordon Smyth.
#	Created 18 Mar 2013. Last modified 13 Apr 2026.
{
#	Try to find gene lengths
#	If column name containing gene lengths isn't specified,
#	then will try "Length" or "length" or any column name containing "length"
	if(is.character(gene.length)) {
		gene.length <- y$genes[[gene.length[1]]]
		if(is.null(gene.length)) stop("gene.length column not found")
	} else {
		if(is.null(gene.length)) gene.length <- y$genes$Length
		if(is.null(gene.length)) gene.length <- y$genes$length
		if(is.null(gene.length)) gene.length <- y$genes$AveLength
		if(is.null(gene.length)) {
			j <- grep("length",tolower(names(y$genes)))
			if(length(j)==1)
				gene.length <- y$genes[,j]
			else
				stop("Gene lengths not found")
		}
	}

	lib.size <- y$samples$lib.size
	if(hasName(y,"offset")){
		if( min(y$offset) > max(log(lib.size)) || min(log(lib.size)) > max(y$offset) ) warning("Offset may not reflect library sizes. Scaling offset may be required.")
		lib.size <- NULL
	} else {
		if(normalized.lib.sizes) lib.size <- lib.size*y$samples$norm.factors
	}

	rpkm.default(y=y$counts, gene.length=gene.length, lib.size=lib.size, offset=y[["offset"]], offset.prior=y[["offset.prior"]], log=log, prior.count=prior.count, ...)
}

rpkm.SummarizedExperiment <- function(y, gene.length=NULL, normalized.lib.sizes=TRUE, log=FALSE, prior.count=2, ...)
#	RPKM for a SummarizedExperiment.
#	Created 03 Apr 2020.  Last modified 03 Apr 2020.
{
	y <- SE2DGEList(y)
	rpkm.DGEList(y, gene.length=gene.length, normalized.lib.sizes=normalized.lib.sizes, log=log, prior.count=prior.count, ...)
}

rpkm.DGELRT <- rpkm.DGEGLM <- rpkm.MArrayLM <- function(y, gene.length, log=FALSE, shrunk = TRUE, ...)
#	Fitted RPKM from a limma or edgeR fitted model object.
#	Gordon Smyth
#	Created 1 Apr 2020. Last modified 9 Oct 2025.
{
#	Try to find gene lengths
#	If column name containing gene lengths isn't specified,
#	then will try "Length" or "length" or any column name containing "length"
	if(is.character(gene.length)) {
		gene.length <- y$genes[[gene.length[1]]]
		if(is.null(gene.length)) stop("gene.length column not found")
	} else {
		if(is.null(gene.length)) gene.length <- y$genes$Length
		if(is.null(gene.length)) gene.length <- y$genes$length
		if(is.null(gene.length)) {
			j <- grep("length",tolower(names(y$genes)))
			if(length(j)==1)
				gene.length <- y$genes[,j]
			else
				stop("Gene lengths not found")
		}
	}

	y <- cpm(y, gene.length=gene.length, log=log, shrunk=shrunk, ...)
	gene.length.kb <- gene.length/1000
	if(log)
		y-log2(gene.length.kb)
	else
		y/gene.length.kb
}

rpkm.default <- function(y, gene.length, lib.size=NULL, offset=NULL, offset.prior=NULL, log=FALSE, prior.count=2, ...)
#	Reads per kilobase of gene length per million reads of sequencing (RPKM)
#	Gordon Smyth
#	Created 1 Nov 2012. Last modified 12 Apr 2026.
{
	y <- cpm.default(y=y, lib.size=lib.size, offset=offset, offset.prior=offset.prior, log=log, prior.count=prior.count, ...)
	gene.length.kb <- gene.length/1000
	if(log)
		y-log2(gene.length.kb)
	else
		y/gene.length.kb
}

