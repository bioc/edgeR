mglmLevenberg <- function(y, design, dispersion=0, offset=0, weights=NULL, coef.start=NULL, start.method="null", maxit=200, tol=1e-06)
#	Fit genewise negative binomial glms with log-link
#	using Levenberg damping to ensure convergence

#	R version by Gordon Smyth and Yunshun Chen
#	C++ version by Aaron Lun
#	C version by Lizhong Chen
#	Created 3 March 2011.  Last modified 16 Sep 2024
{
#	Check arguments
	y <- as.matrix(y)
	if(!is.numeric(y)) stop("y is non-numeric")
	nlibs <- ncol(y)
	ngenes <- nrow(y)
	if(identical(nlibs,0L) || identical(ngenes,0L)) stop("no data")

#	Check for NA, negative or infinite counts
	m <- min(y)
	if(is.na(m)) stop("NA counts not allowed")
	if(m < 0) stop("Negative counts not allowed")
	if(!is.finite(max(y))) stop("Infinite counts not allowed")

#	Checking the design matrix
	design <- as.matrix(design)
	if (!is.double(design)) storage.mode(design) <- "double"
	if (!all(is.finite(design))) stop("all entries of design matrix must be finite and non-missing")

#	Checking dispersions, offsets and weights
	offset <- .compressOffsets(y, offset=offset)
    dispersion <- .compressDispersions(y, dispersion)
	weights <- .compressWeights(y, weights)

#	Checking empty design matrix
	if(identical(ncol(design),0L)){
		output <- list()
		output$coefficients  <- matrix(0,ngenes,0)
		output$fitted.values <- exp(as.matrix(offset))
		output$deviance      <- nbinomDeviance(y=y,mean=output$fitted.values,dispersion=dispersion,weights=weights)
		output$iter          <- rep_len(0L,ngenes)
		output$failed        <- rep_len(FALSE,ngenes) 

		colnames(output$coefficients)  <- colnames(design)
		rownames(output$coefficients)  <- rownames(y)
		dimnames(output$fitted.values) <- dimnames(y)
		
		return(output)
	}

#	Initializing values for the coefficients at reasonable best guess with linear models.
	if(is.null(coef.start)) {
		start.method <- match.arg(start.method, c("null","y"))
		beta <- .Call(.cxx_get_levenberg_start, y, offset, dispersion, weights, design, start.method=="null")
	} else {
		beta <- as.matrix(coef.start)
	}

# 	Checking arguments and calling the C++ method.
	if (!is.double(beta)) storage.mode(beta) <- "double"
	output <- .Call(.cxx_fit_levenberg, y, offset, dispersion, weights, design, beta, tol, maxit)

#	Naming the output and returning it.  
	colnames(output$coefficients) <- colnames(design)
	rownames(output$coefficients) <- rownames(y)
	output
}
