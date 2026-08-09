cpm <- function(y, ...)
UseMethod("cpm")

cpm.DGEList <- function(y, normalized.lib.sizes=TRUE, log=FALSE, prior.count=2, ...)
#	Counts per million for a DGEList
#	Davis McCarthy and Gordon Smyth.
#	Created 20 June 2011. Last modified 9 Aug 2026.
{
	if(hasName(y,"offset"))
		return(cpm(y$counts, offset=y[["offset"]], log=log, prior.count=prior.count))

	lib.size <- y$samples$lib.size
	if(normalized.lib.sizes) lib.size <- lib.size*y$samples$norm.factors
	cpm(y$counts, lib.size=lib.size, offset.prior=y[["offset.prior"]], log=log, prior.count=prior.count)
}

cpm.SummarizedExperiment <- function(y, normalized.lib.sizes=TRUE, log=FALSE, prior.count=2, ...)
#	Counts per million for a SummarizedExperiment
#	Created 03 April 2020.  Last modified 1 June 2020.
{
	y <- SE2DGEList(y)
	cpm.DGEList(y, normalized.lib.sizes=normalized.lib.sizes, log=log, prior.count=prior.count, ...)
}

cpm.DGELRT <- cpm.DGEGLM <- function(y, log=FALSE, shrunk=TRUE, ...)
#	Fitted counts per million from a fitted model object.
#	Created 19 April 2020.  Last modified 30 Sep 2025.
{
	if(shrunk) {
		eta <- y$coefficients %*% t(y$design)
	} else {
		eta <- y$unshrunk.coefficients %*% t(y$design)
	}

	if(log) {
		(eta + log(1e6)) / log(2)
	} else {
		exp(eta + log(1e6))
	}
}

cpm.MArrayLM <- function(y, log=FALSE, ...)
#	Fitted counts per million from a limma fitted model object.
#	Created 9 Oct 2025.  Last modified 9 Oct 2025.
{
	eta <- fitted(y)
	if(log) {
		eta
	} else {
		2^eta
	}
}

cpm.default <- function(y, lib.size=NULL, offset=NULL, offset.prior=NULL, log=FALSE, prior.count=2, ...)
#	Counts per million for a matrix
#	Davis McCarthy, Yunshun Chen, Gordon Smyth.
#   C++ version by Aaron Lun. C version by Lizhong Chen.
#	Created 20 June 2011. Last modified 9 Aug 2026.
{
#	Coerce to matrix
	y <- as.matrix(y)

#	Check for zero length
	if(!length(y)) return(y)

#	Check y entries
	ymin <- min(y)
	if(is.na(ymin)) stop("NA counts not allowed")
	if(ymin < 0) stop("Negative counts not allowed")

#	Offset matrix takes precedence over lib.size and offset.prior. Otherwise, lib.sizes default to column sums.
#	offset.prior adds to log(lib.size) while offset replaces them.
	if(is.null(offset)) {
		if(is.null(lib.size)) lib.size <- colSums(y)
		if(!is.null(offset.prior)) {
			if(!identical(dim(y),dim(offset.prior))) stop("counts and offset.prior must have equal dimensions.")
			lib.size <- exp(matrix(log(lib.size),nrow(y),ncol(y),byrow=TRUE) + offset.prior)
		}
	} else {
#		Offset can be a matrix or a row vector.
		if(is.matrix(offset)) {
			if(any(dim(offset)!=dim(y))) stop("counts and offset must have equal dimensions")
		} else {
			if(length(offset)!=ncol(y)) stop("if offset is a vector, its length must be the number of samples")
		}
		if(!is.null(lib.size)) warning("lib.size is ignored in the presence of offset")
		if(!is.null(offset.prior)) warning("offset.prior is ignored in the presence of offset")
		lib.size <- exp(offset)
	}

#	Ensure lib.size is double precision numeric for compatibility with C
	if(!is.double(lib.size)) {
		if(!is.numeric(lib.size)) stop("lib.size must be numeric")
		storage.mode(lib.size) <- "double"
	}

#	Check lib.sizes are not missing or zero
	minlibsize <- min(lib.size)
	if(is.na(minlibsize)) stop("NA library sizes not allowed")
	if(minlibsize <= 0) stop("Library sizes should be greater than zero")

#	Convert to compressed matrix, allowing for lib.size as either a matrix or a row vector
	lib.size <- makeCompressedMatrix(lib.size, dim(y), byrow=TRUE)

#	Calculating in C for max efficiency
	if(log) {
		prior.count <- .compressPrior(y, prior.count)
		out <- .Call(.cxx_calculate_cpm_log, y, lib.size, prior.count)
	} else {
		out <- .Call(.cxx_calculate_cpm_raw, y, lib.size)
	}

	out
}
