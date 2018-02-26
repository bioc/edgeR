cpmByGroup <- function(y, ...)
UseMethod("cpmByGroup")

cpmByGroup.DGEList <- function(y, group=NULL, dispersion=NULL, ...)
#	Counts per million averaged by group
#	Gordon Smyth
#	Created 10 July 2017. Last modified 5 Oct 2017.
{
	if(is.null(group)) group <- y$samples$group
	group <- as.factor(group)

	if(is.null(dispersion)) dispersion <- getDispersion(y)
	if(is.null(dispersion)) dispersion <- 0.05
	offset <- getOffset(y)

	fit <- mglmOneWay(y,group=group,dispersion=dispersion,offset=offset,weights=y$weights)
	exp(fit$coefficients) * 1e6
}

cpmByGroup.default <- function(y, group=NULL, dispersion=0.05, offset=NULL, weights=NULL, ...)
#	Counts per million averaged by group
#	Gordon Smyth
#	Created 10 July 2017. Last modified 26 Feb 2018.
{
	y <- as.matrix(y)

	if(is.null(group)) {
		group <- factor(rep_len(1,ncol(y)))
		levels(group) <- "AveCPM"
	}

	if(is.null(offset)) offset <- log(colSums(y))

	fit <- mglmOneWay(y,group=group,dispersion=dispersion,offset=offset,weights=weights)
	exp(fit$coefficients) * 1e6
}
