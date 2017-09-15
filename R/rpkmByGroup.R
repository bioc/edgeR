rpkmByGroup <- function(y, ...)
UseMethod("rpkmByGroup")

rpkmByGroup.DGEList <- function(y, group=NULL, gene.length=NULL, log=FALSE, prior.count=0.25, ...)
#	Counts per million averaged by group
#	Gordon Smyth
#	Created 10 July 2017
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

	if(!is.null(group)) y$samples$group <- group
	if(!log) prior.count <- 0
	fit <- glmFit(y,prior.count=prior.count,...)

	RPKM <- fit$coefficients - log(gene.length) + log(1e9)
	if(log) {
		RPKM / log(2)
	} else {
		exp(RPKM)
	}
}

rpkmByGroup.default <- function(y, group=NULL, gene.length=NULL, dispersion=0.05, log=FALSE, prior.count=0.25, ...)
#	Counts per million averaged by group
#	Gordon Smyth
#	Created 10 July 2017
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

	if(!log) prior.count <- 0
	if(is.null(group)) {
		fit <- glmFit(y,dispersion=dispersion,prior.count=prior.count,...)
	} else {
		group <- as.factor(group)
		design <- model.matrix(~0+group)
		colnames(design) <- levels(group)
		fit <- glmFit(y,design=design,dispersion=dispersion,prior.count=prior.count,...)
	}

	logRPKM <- fit$coefficients - log(gene.length) + log(1e9)
	if(log) {
		logRPKM / log(2)
	} else {
		exp(logRPKM)
	}
}
