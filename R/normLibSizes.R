# This file implements a renaming of calcNormFactors() to normLibSizes() on 7 Nov 2022

normLibSizes <- function(object, ...)
UseMethod("normLibSizes")

normLibSizes.DGEList <- function(object, method=c("TMM","TMMwsp","RLE","upperquartile","none"), refColumn=NULL, logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10, p=0.75, renormalize.offset=FALSE, ...)
#	Library size normalization for read count matrices.
#	Method for DGEList objects.
#	Created 2 Oct 2014.
#	Support for offset.prior added 1 Nov 2026.
#	Last modified 1 Nov 2026.
{
	if(hasName(object,"offset")) {
		if(hasName(object,"offset.prior")) stop("Object contains both offset and offset.prior matrices.")
		if(renormalize.offset) {
			message("Converting offset matrix to offset.prior with normalized library sizes")
			object$offset.prior <- object$offset - matrix(log(object$samples$lib.size),nrow(object),ncol(object),byrow=TRUE)
		} else {
			message("Data object contains offset matrix so skipping library size normalization")
			return(object)
		}
	}
	object$samples$norm.factors <- normLibSizes(
		object=object$counts,
		lib.size=object$samples$lib.size,
		offset.prior=object$offset.prior,
		method=method,
		refColumn=refColumn,
		logratioTrim=logratioTrim,
		sumTrim=sumTrim,
		doWeighting=doWeighting,
		Acutoff=Acutoff,
		p=p)
	object
}

normLibSizes.SummarizedExperiment <- function(object, method=c("TMM","TMMwsp","RLE","upperquartile","none"), refColumn=NULL, logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10, p=0.75, ...)
#	Library size normalization for read count matrices.
#	Method for SummarizedExperiment objects.
#	Created 19 Mar 2020. Last modified 7 Nov 2022.
{
	object <- SE2DGEList(object)
	object$samples$norm.factors <- normLibSizes(object=object$counts, lib.size=object$samples$lib.size, method=method, refColumn=refColumn, logratioTrim=logratioTrim, sumTrim=sumTrim, doWeighting=doWeighting, Acutoff=Acutoff, p=p)
	object
}

normLibSizes.default <- function(object, lib.size=NULL, offset.prior=NULL, method=c("TMM","TMMwsp","RLE","upperquartile","none"), refColumn=NULL, logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10, p=0.75, ...)
#	Library size normalization for read count matrices.
#	Method for atomic count matrices.
#	Mark Robinson, Gordon Smyth and edgeR team
#	Created 22 Oct 2009.
#	Support for offset.prior added 1 Nov 2026.
#	Last modified 1 Nov 2026.
{
#	Check object
	x <- as.matrix(object)
	minx <- min(x)
	if(is.na(minx)) stop("NA counts not allowed")
	if(minx < 0) stop("Negative counts not allowed")
	ngenes <- nrow(x)
	nsamples <- ncol(x)

#	Check lib.size
	if(is.null(lib.size)) {
		lib.size <- colSums(x)
	} else {
		if(anyNA(lib.size)) stop("NA lib.sizes not permitted")
		if(length(lib.size) != nsamples) {
			if(length(lib.size) > 1L) warning("normLibSizes: length(lib.size) doesn't match number of samples",call.=FALSE)
			lib.size <- rep_len(lib.size,nsamples)
		}
	}

#	Combine offset.prior with library sizes
	if(!is.null(offset.prior)) {
		if(!identical(dim(object),dim(offset.prior))) stop("data object and offset.prior are not of same dimensions")
		lib.size.matrix <- exp( offset.prior + matrix(log(lib.size),ngenes,nsamples,byrow=TRUE) )
	}

#	Check method
#	Backward compatability with previous name
	if(length(method)==1L && method=="TMMwzp") {
		method <- "TMMwsp"
		message("TMMwzp has been renamed to TMMwsp")
	}
	method <- match.arg(method)
	if(!is.null(offset.prior) && method!="TMM") message("ignoring offset.prior for ",method," method")

#	Remove all zero rows
	allzero <- .rowSums(x>0, nrow(x), nsamples) == 0L
	if(any(allzero)) x <- x[!allzero,,drop=FALSE]

#	Degenerate cases
	if(nrow(x)==0 || nsamples==1) method="none"

#	Calculate factors
	f <- switch(method,
		TMM = {
			if( is.null(refColumn) ) {
				f75 <- rep_len(1,nsamples)
				if(is.null(offset.prior)) {
					for(i in 1:nsamples) f75[i] <- quantile(x[,i], probs=0.75) / lib.size[i]
				} else {
					for(i in 1:nsamples) f75[i] <- quantile(x[,i]/lib.size.matrix[,i], probs=0.75)
				}
				if(median(f75) < 1e-20) {
					refColumn <- which.max(colSums(sqrt(x)))
				} else {
					refColumn <- which.min(abs(f75-mean(f75)))
				}
			}
			f <- rep_len(NA_real_,nsamples)
			Yref <- x[,refColumn]
			if(is.null(offset.prior)) {
				Nref <- lib.size[refColumn]
				for(i in 1:nsamples) f[i] <- .calcFactorTMM(obs=x[,i], ref=Yref, libsize.obs=lib.size[i], libsize.ref=Nref, logratioTrim=logratioTrim, sumTrim=sumTrim, doWeighting=doWeighting, Acutoff=Acutoff)
			} else {
				Nref <- lib.size.matrix[,refColumn]
				for(i in 1:nsamples) f[i] <- .calcFactorTMM(obs=x[,i], ref=Yref, libsize.obs=lib.size.matrix[,i], libsize.ref=Nref, logratioTrim=logratioTrim, sumTrim=sumTrim, doWeighting=doWeighting, Acutoff=Acutoff)
			}
			f
		},
		TMMwsp = {
			if( is.null(refColumn) ) refColumn <- which.max(colSums(sqrt(x)))
			f <- rep_len(NA_real_,nsamples)
			for(i in 1:nsamples)
				f[i] <- .calcFactorTMMwsp(obs=x[,i], ref=x[,refColumn], libsize.obs=lib.size[i], libsize.ref=lib.size[refColumn], logratioTrim=logratioTrim, sumTrim=sumTrim, doWeighting=doWeighting, Acutoff=Acutoff)
			f
		},
		RLE = .calcFactorRLE(x)/lib.size,
		upperquartile = .calcFactorQuantile(x,lib.size,p=p),
		none = rep_len(1,nsamples)
	)

#	Factors should multiple to one
	f <- f/exp(mean(log(f)))

#	Output
	names(f) <- colnames(x)
	f
}
