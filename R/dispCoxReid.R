dispCoxReid <- function(y, design=NULL, offset=NULL, weights=NULL, AveLogCPM=NULL, interval=c(0,4), tol=1e-5, min.row.sum=5, subset=10000, nthreads=1L)
#	Cox-Reid APL estimator of common dispersion
#	Gordon Smyth, Davis McCarthy
#	26 Jan 2011.  Last modified 9 Dec 2013.
{
#	Check y
	y <- as.matrix(y)

#	Check design
	if(is.null(design)) {
		design <- matrix(1,ncol(y),1)
		rownames(design) <- colnames(y)
		colnames(design) <- "Intercept"
	} else {
		design <- as.matrix(design)
	}

#	Check offset
	if(is.null(offset)) offset <- log(colSums(y))
	offset <- expandAsMatrix(offset,dim(y))
	if(min(interval)<0) stop("please give a non-negative interval for the dispersion")

#	Apply min row count
	small.row.sum <- rowSums(y)<min.row.sum
	if(any(small.row.sum)) {
		y <- y[!small.row.sum,,drop=FALSE]
		offset <- offset[!small.row.sum,,drop=FALSE]
		weights <- weights[!small.row.sum,,drop=FALSE]
		if(!is.null(AveLogCPM)) AveLogCPM <- AveLogCPM[!small.row.sum]
	}
	if(nrow(y)<1) stop("no data rows with required number of counts")

#	Subsetting
	if(!is.null(subset) && subset<=nrow(y)/2) {
		if(is.null(AveLogCPM)) AveLogCPM <- aveLogCPM(y,offset=offset,weights=weights)
		i <- systematicSubset(subset,AveLogCPM)
		y <- y[i,,drop=FALSE]
		offset <- offset[i,,drop=FALSE]
		weights <- weights[i,,drop=FALSE]
	}

#	Validate the design and compress inputs once (these were previously repeated
#	on every optimize() evaluation inside adjustedProfileLik/glmFit).
	ne <- nonEstimable(design)
	if(!is.null(ne)) stop(paste("Design matrix not of full rank. The following coefficients not estimable:\n", paste(ne, collapse=" ")))
	offset  <- .compressOffsets(y, offset=offset)
	weights <- .compressWeights(y, weights)

#	Maximize the summed Cox-Reid adjusted profile likelihood over the dispersion
#	in a single C call: a 1-D Brent search (in par = disp^0.25 space) whose
#	objective fits all genes and sums the APL, warm-starting each fit from the
#	previous dispersion. Returns the optimal dispersion (par^4) directly.
	.Call(.cxx_coxreid_disp, y, offset, weights, design,
	      min(interval)^0.25, max(interval)^0.25, tol, TRUE,
	      nthreads)
}
