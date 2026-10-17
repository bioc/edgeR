#  sumTechReps.R

sumTechReps <- function(x,ID=colnames(x),...) UseMethod("sumTechReps")
#  24 Sept 2010

sumTechReps.default <- function(x,ID=colnames(x),...)
#	Sum over replicate columns, for matrices
#	Yifang Hu and Gordon Smyth
#	Created 14 Mar 2014
{
	if(is.null(x)) return(NULL)
	x <- as.matrix(x)
	if(is.null(ID)) stop("No sample IDs")
	t(rowsum(t(x),group=ID,reorder=FALSE,na.rm=FALSE))
}

sumTechReps.DGEList <- function(x,ID=colnames(x),...)
#	Sum over replicate columns, for DGELists
#	Gordon Smyth and Yifang Hu
#	Created 14 Mar 2014. Last modified 17 Oct 2026.
{
#	If no tech reps, just return object
	d <- duplicated(ID)
	if(!any(d)) return(x)

#	Remove any existing dispersion estimates
	x$common.dispersion <- x$trended.dispersion <- x$tagwise.dispersion <- NULL

#	Weights are not currently handled
	if(!is.null(x$weights)) {
		warning("Removing weights, which are not currently handled by sumTechReps")
		x$weights <- NULL
	}

#	Collapse object to unique sample IDs
	y <- x[,!d]

#	Count number of tech reps
	sumReps <- function(z,ID) rowsum(z,group=ID,reorder=FALSE,na.rm=FALSE)
	n <- drop(sumReps(rep_len(1L,ncol(x)),ID))

#	Sum counts
	y$counts <- t(sumReps(t(x$counts),ID))

#	Sum library sizes
	y$samples$lib.size <- drop(sumReps(x$samples$lib.size,ID))

#	Geometric mean of normalization factors
	y$samples$norm.factors <- exp(drop(sumReps(log(x$samples$norm.factors),ID))/n)

#	Sum offsets on unlogged scale
	if(hasName(x,"offset")) {
		y$offset <- log(t(sumReps(exp(t(x[["offset"]])),ID)))
		if(hasName(x,"offset.prior")) warning("Object has both offset and offset.prior. Latter will be ignored.")
	}

#	Average offset.prior matrix
	if(hasName(x,"offset.prior")) {
		y$offset.prior <- t(rowsum(t(x[["offset.prior"]]),group=ID,reorder=FALSE,na.rm=FALSE)/n)
	}

#	New sample names
	rownames(y$samples) <- colnames(y$counts)

	y
}

sumTechReps.SummarizedExperiment <- function(x, ID, ...)
#	Sum over technical replicates, for SummarizedExperiment objects
#	Created 03 Apr 2020.  Last modified 03 Apr 2020.
{
	x <- SE2DGEList(x)
	sumTechReps.DGEList(x, ID=ID, ...)
}
