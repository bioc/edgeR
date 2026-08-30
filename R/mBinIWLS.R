mbinIWLS <- function (y, z, design, offset = NULL, weights = NULL, maxit = 100, tol = 1e-06, nthreads=1L)
#	Fit binomial logistic models using iterated weighted least square
#	Created by Lizhong Chen
#	11 April 2025.  Last modified 15 April 2025.
{
#	check counts
	y <- as.matrix(y)
	z <- as.matrix(z)
	
	ngenes <- nrow(y)
	nlibs  <- ncol(y)

#	check design
	design <- as.matrix(design)
	if(nrow(design) != nlibs) stop("nrow(design) disagrees with ncol(y)")
	if (!all(is.finite(design))) stop("all entries of design matrix must be finite and non-missing")

#	check weights
	coverage <- y + z
	if(is.null(weights)){
		weights <- coverage    
	} else {
		weights <- .compressWeights(y, weights)
		weights <- as.matrix(weights) * coverage
	}
	
#	calculate proportion
	prop <- y / pmax(coverage,1)
	
#	check offset
	if(is.null(offset)) {
		offset <- makeCompressedMatrix(0,dim(y))
	} else {
		offset <- .compressOffsets(y, offset)
	}

	out <- .Call(.cxx_bin_fit_iwls, prop, offset, weights, design, maxit, tol, nthreads)

	dimnames(out$coefficients) <- list(rownames(y),colnames(design))

	out
}
