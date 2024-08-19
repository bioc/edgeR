loessByCol <- function(y, x=NULL, span=0.5)
# Fit a lowess curve of degree 0 by least squares to each column of a matrix.
#
# R version called locallyWeightedMean() using simpleLoess() by Yunshun Chen, 08 May 2012.
# Converted to C++ and renamed to loessByCol() by Aaron Lun, 26 June 2012.
# Pure C version by Lizhong Chen, 11 May 2024.
{
#	Check y
	y <- as.matrix(y)
	ntags <- nrow(y)

#	Check x
	if(is.null(x)) x <- 1:ntags
	if(!identical(length(x),ntags)) stop("Length of `x` must match row dimension of 'y'")
	names(x) <- rownames(y)

#	If span window is less than one observation, return y without smoothing
	nspan <- min(floor(span*ntags), ntags)
	if(nspan<=1) {
		fitted <- list(fitted.values=y,leverages=rep(1,ntags))
		names(fitted$leverages) <- rownames(y)
		return(fitted)
	}

#	Call C code
	.Call(.cxx_loess_by_col, x, y, nspan)
}
