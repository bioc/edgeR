locfitByCol <- function(y, x=NULL, weights=1, span=0.5, degree=0, nthreads=1L)
#	
#	Gordon Smyth
#	Created 20 Aug 2012. Last modified 30 Oct 2026.
{
	y <- as.matrix(y)
	ntags <- nrow(y)
	weights <- rep_len(weights,ntags)
	if(is.null(x)) x <- 1:ntags
	if(length(degree)!=1L || !(degree %in% 0:3)) stop("degree must be 0, 1, 2 or 3")
	if(span*ntags<2 || ntags<=1 || ncol(y)==0L) return(y)
	if(!all(is.finite(y)) || !all(is.finite(x)) || !all(is.finite(weights))) stop("x, y and weights must be finite")
#	fitted(locfit(y[,j]~x, weights=weights, alpha=span, deg=degree)) for each column, in C (locfit_lite.c)
	.Call(.cxx_locfit_by_col, y, x, weights, span, degree, nthreads)
}
