#  FIT GENERALIZED LINEAR MODELS

filterByExpr <- function(y, ...)
UseMethod("filterByExpr")

filterByExpr.DGEList <- function(y, design=NULL, group=NULL, lib.size=NULL, ...)
#	Filter low expressed genes given count matrix
#	Gordon Smyth
#	Created 13 Nov 2017. Last revised 9 Aug 2026.
{
#	Order of precedence:
#	1. group as argument
#	2. design as argument
#	3. y$design
#	4. y$samples$group
	if(is.null(design) && is.null(group)) {
		design <- y$design
		if(is.null(design)) {
			group <- y$samples$group
			if(length(levels(group))==1L) message("Treating all samples as belonging to one group.")
		}
	}
	if(is.null(lib.size)) lib.size <- y$samples$lib.size * y$samples$norm.factors
	filterByExpr.default(y$counts, design=design, group=group, lib.size=lib.size, ...)
}

filterByExpr.SummarizedExperiment <- function(y, design=NULL, group=NULL, lib.size=NULL, ...)
#	Created 19 Mar 2020. Last revised 19 Mar 2020.
{
	y <- SE2DGEList(y)
	filterByExpr.DGEList(y, design=design, group=group, lib.size=lib.size, ...)
}

filterByExpr.default <- function(y, design=NULL, group=NULL, min.sample.size=NULL, lib.size=NULL, min.count=10, min.total.count=15, large.n=10, min.prop=0.7, max.n.samples=40, preset=NULL, ...)
#	Filter low expressed genes given count matrix.
#	Computes TRUE/FALSE index vector indicating which rows to keep.
#	Gordon Smyth
#	Created 13 Nov 2017.
#	Arguments min.n.samples, max.n.samples, and preset added 9 Oct 2026.
#	Last modified 18 Oct 2026.
{
	y <- as.matrix(y)
	if(mode(y) != "numeric") stop("y is not a numeric matrix")
	if(is.null(lib.size)) lib.size <- colSums(y)

#	Preset options.
#	The "stringent" option matches the historical default setting for this function from 2020 to now.
#	The "lenient" option does minimal filtering, assuming that
#	functions like glmQLFit() and voomLmFit() are now robust to small counts.
	if(!is.null(preset)) {
		preset <- match.arg(preset, c("stringent","moderate","lenient"))
		if(!all(missing(min.count),missing(min.total.count),missing(large.n),missing(min.prop),missing(max.n.samples),missing(max.n.samples)))
			message("Using ",preset," settings. Ignoring arguments set individually.")
		switch(preset,
		stringent = {
			min.count <- 10
			min.total.count <- 15
			large.n <- 10
			min.prop <- 0.7
			max.n.samples <- 40
		},
		moderate = {
			min.count <- 6
			min.total.count <- 10
			large.n <- 4
			min.prop <- 0.5
			max.n.samples <- 20
		},
		lenient = {
			min.count <- 1
			min.total.count <- 6
			large.n <- 2
			min.prop <- 0.5
			max.n.samples <- 10
		})
	}

#	Minimum effect sample sample size for any of the coefficients
	if(is.null(min.sample.size)) {
		if(is.null(group)) {
			if(is.null(design)) {
				message("No group or design set. Assuming all samples belong to one group.")
				MinSampleSize <- ncol(y)
			} else {
				h <- hat(design)
				MinSampleSize <- 1/max(h)
			}
		} else {
			group <- as.factor(group)
			n <- tabulate(group)
			MinSampleSize <- min(n[n > 0L])
		}
	} else {
		MinSampleSize <- min.sample.size
	}

#	For the first few samples, require genes to be expressed in every one.
#	As the group sample size increases, require genes to be expressed in just a proportion of the additional samples.
	if(MinSampleSize > large.n) MinSampleSize <- large.n + (MinSampleSize-large.n)*min.prop

#	When all the group sizes are large, keep any gene expressed in max.n.samples samples
	MinSampleSize <- min(MinSampleSize, max.n.samples)

#	CPM cutoff
	MedianLibSize <- median(lib.size)
	CPM.Cutoff <- min.count/MedianLibSize*1e6
	CPM <- cpm(y,lib.size=lib.size)
	tol <- 1e-14
	keep.CPM <- rowSums(CPM >= CPM.Cutoff) >= (MinSampleSize - tol)

#	Total count cutoff
	keep.TotalCount <- (rowSums(y) >= min.total.count - tol)

	keep.CPM & keep.TotalCount
}
