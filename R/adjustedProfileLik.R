adjustedProfileLik <- function(dispersion, y, design, offset, weights=NULL, adjust=TRUE, start=NULL, get.coef=FALSE, nthreads=1L)
#	Tagwise Cox-Reid adjusted profile log-likelihoods for the dispersion.
#	dispersion can be a scalar or a tagwise vector.
#	Computationally, dispersion can also be a matrix, but the apl is still computed tagwise.
#	y is a matrix: rows are genes/tags/transcripts, columns are samples/libraries.
#	offset is a matrix of the same dimensions as y.

#	The weights argument was added by Xiaobei Zhou 20 March 2013,
#	but the log NB probabilities were incorrectly multiplied by the weights.
#	This is fixed 1 March 2018 with a more rigorous interpretation of weights
#	in terms of averages.

#	Yunshun Chen, Gordon Smyth, Aaron Lun
#	Created June 2010. Last modified 22 May 2020.
{
#	Checking counts
	if (!is.numeric(y)) stop("counts must be numeric")
	y <- as.matrix(y)

#	Checking design (full rank); this guard was previously provided by glmFit
	design <- as.matrix(design)
	ne <- nonEstimable(design)
	if(!is.null(ne)) stop(paste("Design matrix not of full rank. The following coefficients not estimable:\n", paste(ne, collapse=" ")))

#	Checking offsets
	offset <- .compressOffsets(y, offset=offset)

#	Checking dispersion
	dispersion <- .compressDispersions(y, dispersion)

#	Checking weights
	weights <- .compressWeights(y, weights)
	  
#	Fit tagwise GLMs and compute the adjusted profile likelihood in one C call.
#	fit_glm_mat produces the fitted values internally and feeds them straight to
#	the APL computation, so 'mu' never round-trips through R. The 250L/1e-6 match
#	the maxit/tol that glmFit.default passes to .cxx_fit_glm.
	fit <- .Call(.cxx_fit_apl, y, offset, dispersion, weights, design, 250L, 1e-6, start, adjust, nthreads)

#	Deciding what to return.
	if (get.coef) {
		beta <- fit$coefficients
		dimnames(beta) <- list(rownames(y), colnames(design))
		return(list(apl=fit$apl, beta=beta))
	} else {
		return(fit$apl)
	}
}

