cpm <- function(y, ...)
UseMethod("cpm")

cpm.DGEList <- function(y, normalized.lib.sizes=TRUE, log=FALSE, prior.count=2, ...)
#	Counts per million for a DGEList
#	Davis McCarthy and Gordon Smyth.
#	Created 20 June 2011. Last modified 12 Apr 2026.
{
	lib.size <- y$samples$lib.size
	if(normalized.lib.sizes) lib.size <- lib.size*y$samples$norm.factors
	cpm.default(y$counts, lib.size=lib.size, offset=y[["offset"]], offset.prior=y[["offset.prior"]], log=log, prior.count=prior.count)
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
	# check paired counts
	if(is.null(y$dispersion)) stop("cpm does not work for paired counts")

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
#	Davis McCarthy and Gordon Smyth. C++ version by Aaron Lun. C version by Lizhong Chen.
#	Created 20 June 2011. Last modified 12 Apr 2026.
{
#	Check y
	ymin <- min(y)
	if(is.na(ymin)) stop("NA counts not allowed")
	if(ymin < 0) stop("Negative counts not allowed")
	y <- as.matrix(y)
	if(any(dim(y)==0L)) {
		return(y)
	}
#	Check lib.size
	if(is.null(lib.size)) lib.size <- colSums(y)

#	Check lib.sizes are not missing or zero
	minlibsize <- min(lib.size)
	if(is.na(minlibsize)) stop("NA library sizes not allowed")
	if(minlibsize <= 0) stop("Library sizes should be greater than zero")

#	Ensure lib.size is double precision numeric for compatibility with C
	if(!is.double(lib.size)) {
		if(!is.numeric(lib.size)) stop("lib.size must be numeric")
		storage.mode(lib.size) <- "double"
	}

#	Combine library sizes with prior offset matrix (if present)
	if(!is.null(offset)) {
		if(is.null(offset.prior)) {
			if(!identical(dim(y),dim(offset))) stop("y and offset must have equal dimensions.")
			offset.prior <- offset - rowMeans(offset)
		} else {
			message("Ignoring offset in favor of offset.prior. Should not set both.")
			offset <- NULL
		}
	}
	if(!is.null(offset.prior)) {
		if(!identical(dim(y),dim(offset.prior))) stop("y and offset.prior must have equal dimensions.")
		lib.size <- exp(matrix(log(lib.size),nrow(y),ncol(y),byrow=TRUE) + offset.prior)
	}

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
