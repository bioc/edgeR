rpkm <- function(y, ...)
UseMethod("rpkm")

rpkm.DGEList <- function(y, gene.length=NULL, normalized.lib.sizes=TRUE, log=FALSE, prior.count=0.25, ...)
#	Reads per kilobase of gene length per million reads of sequencing (RPKM)
#	Gordon Smyth.
#	Created 18 March 2013. Last modified 1 November 2012
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

	lib.size <- y$samples$lib.size
	if(normalized.lib.sizes) lib.size <- lib.size*y$samples$norm.factors

	rpkm.default(y=y$counts,gene.length=gene.length,lib.size=lib.size,log=log,prior.count=prior.count, ...)
}

rpkm.default <- function(y, gene.length, lib.size=NULL, log=FALSE, prior.count=0.25, ...)
#	Reads per kilobase of gene length per million reads of sequencing (RPKM)
#	Gordon Smyth
#	Created 1 November 2012. Last modified 18 March 2014.
{
	y <- cpm.default(y=y,lib.size=lib.size,log=log,prior.count=prior.count, ...)
	gene.length.kb <- gene.length/1000
	if(log)
		y-log2(gene.length.kb)
	else
		y/gene.length.kb
}

