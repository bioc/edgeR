nbinomDeviance <- function(y,mean,dispersion=0,weights=NULL)
#	Negative binomial residual deviance
#	y is a matrix and a deviance is computed for each row
#	A vector y is taken to be a matrix with one row; in this case mean must also be a vector.
#	Original version 23 November 2010.
#	Last modified 21 June 2017.
{
	out <- .compute_nbdeviance(y=y, mean=mean, dispersion=dispersion, weights=weights, dosum=TRUE)
	names(out) <- rownames(y)
	out
}

.compute_nbdeviance <- function(y, mean, dispersion, weights, dosum)
{
#	Check y. May be matrix or vector.
	if(is.matrix(y)) {
		if(!is.matrix(mean)) stop("y is a matrix but mean is not")
	} else {
		n <- length(y)
		y <- matrix(y,1L,n)
		if(is.matrix(mean)) {
			stop("mean is a matrix but y is not")
		} else {
			if(length(mean)==n || length(mean==1L)) {
				mean <- matrix(mean,1L,n)
			}
		}
	}

#	Check mean
	if(!identical(dim(y),dim(mean))) stop("mean should have same dimensions as y")
	if(!is.double(mean)) storage.mode(mean) <- "double"

#	Check dispersion (can be tagwise (rowwise) or observation-wise).
	dispersion <- .compressDispersions(y, dispersion)

#	Check weights.
	weights <- .compressWeights(y, weights)

#	Compute matrix of unit deviances, or residual deviance per gene, depending on 'dosum'.
	.Call(.cxx_compute_nbdev, y, mean, dispersion, weights, as.logical(dosum))
}

nbinomUnitDeviance <- function(y,mean,dispersion=0) 
#	Unit deviance for the nbinom distribution.
#	Created 9 Dec 2013. Last modified 18 Mar 2018.
{
	y[] <- .compute_nbdeviance(y=y, mean=mean, dispersion=dispersion, weights=NULL, dosum=FALSE)
	y
}
