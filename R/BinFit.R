#	Fit Binomial Models
binFit <- function(y, ...)
	UseMethod("binFit")

binFit.PCList <- function(y, design=NULL, ...)
#	Created 15 April 2025.  Last modified 18 April 2025.
{
	if(is.null(design)) {
		design <- y$design
		if(is.null(design)) {
			group <- droplevels(as.factor(y$samples$group))
			if(nlevels(group) > 1L) design <- model.matrix(~y$samples$group)
		}
	}

	if(is.null(y$AveLogCPM)) y$AveLogCPM <- aveLogCPM(y)
	
	fit <- binFit.default(y=y$counts,z=y$counts2,design=design,weights=y$weights, ...)
	
	fit$samples   <- y$samples
	fit$genes     <- y$genes
	fit$AveLogCPM <- y$AveLogCPM
	fit
}

binFit.default <- function(y, z, design=NULL, weights=NULL, offset=NULL, prior.count=1, nthreads=1L, ...)
#	Fit binomial generalized linear model for each feature
#	Lizhong Chen, Gordon Smyth
#	Created 15 April 2025. Last modified 22 April 2025.
{
	#	Check counts
	y <- as.matrix(y)
	z <- as.matrix(z)

	#	Check dim
	if(!identical(dim(y),dim(z))) stop("The dimensions of 'y' and 'z' are not matched")
	
	ntag <- nrow(y)
	nlib <- ncol(y)
	
	#	Check design
	if(is.null(design)) {
		design <- matrix(1,nlib,1)
		rownames(design) <- colnames(y)
		colnames(design) <- "Intercept"
	} else {
		design <- as.matrix(design)
		if(nrow(design) != nlib) stop("nrow(design) disagrees with ncol(y)")
		ne <- nonEstimable(design)
		if(!is.null(ne)) stop(paste("Design matrix not of full rank.  The following coefficients not estimable:\n", paste(ne, collapse = " ")))
	}

	#	coverage, weights and proportion (computed once, before the fit)
	coverage <- y + z
	if(is.null(weights)){
		weights0 <- coverage
	} else {
		weights0 <- as.matrix(weights) * coverage
	}
	prop <- y / pmax(coverage,1)

	#	compressed offset (do not overwrite 'offset'; it is reused by the prior.count recursion below)
	if(is.null(offset)) offset0 <- makeCompressedMatrix(0, dim(y)) else offset0 <- .compressOffsets(y, offset)

	#	Fit the tagwise binomial GLMs in a SINGLE C call.
	#	The oneway-vs-general decision is made inside C, mirroring glmFit / .cxx_fit_glm.
	fit <- .Call(.cxx_bin_fit, prop, offset0, weights0, design, 50L, 1e-8, 100L, 1e-6, nthreads)
	dimnames(fit$coefficients)  <- list(rownames(y), colnames(design))
	dimnames(fit$fitted.values) <- dimnames(y)
	
	#	Prepare output
	fit$counts  <- y
	fit$counts2 <- z

	if(prior.count > 0) {
		fit$unshrunk.coefficients <- fit$coefficients
		y <- y + ceiling(prior.count) / 2
		z <- z + ceiling(prior.count) / 2
		fit$coefficients <- binFit(y,z,design=design,offset=offset,weights=weights,prior.count=0)$coefficients
	}

	fit$df.residual <- rep(nlib-ncol(design),ntag)
	fit$design <- design
	fit$weights <- weights
	fit$nthreads <- nthreads
	new("DGEBIN",fit)
}
