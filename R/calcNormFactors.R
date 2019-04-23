calcNormFactors <- function(object, ...)
UseMethod("calcNormFactors")

calcNormFactors.DGEList <- function(object, method=c("TMMwsp","TMM","RLE","upperquartile","none"), refColumn=NULL, logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10, p=0.75, ...)
#	Scale normalization of RNA-Seq data, for DGEList objects
#	Created 2 October 2014.  Last modified 27 August 2015.
{
	object$samples$norm.factors <- calcNormFactors(object=object$counts, lib.size=object$samples$lib.size, method=method, refColumn=refColumn, logratioTrim=logratioTrim, sumTrim=sumTrim, doWeighting=doWeighting, Acutoff=Acutoff, p=p)
	object
}

calcNormFactors.default <- function(object, lib.size=NULL, method=c("TMMwsp","TMM","RLE","upperquartile","none"), refColumn=NULL, logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10, p=0.75, ...)
#	Scale normalization of RNA-Seq data, for count matrices
#	Mark Robinson.  Edits by Gordon Smyth.
#	Created October 22 October 2009 by Mark Robinson.
#	Last modified 23 Apr 2019.
{
#	Check object
	x <- as.matrix(object)
	if(any(is.na(x))) stop("NA counts not permitted")

#	Check lib.size
	if(is.null(lib.size)) lib.size <- colSums(x)
	if(any(is.na(lib.size))) stop("NA lib.sizes not permitted")

#	Check method
#	Backward compatability with previous name
	if(length(method)==1L && method=="TMMwzp") method="TMMwsp"
	method <- match.arg(method)

#	Remove all zero rows
	allzero <- .rowSums(x>0, nrow(x), ncol(x)) == 0
	if(any(allzero)) x <- x[!allzero,,drop=FALSE]

#	Degenerate cases
	if(nrow(x)==0 || ncol(x)==1) method="none"

#	Calculate factors
	f <- switch(method,
		TMM = {
			f75 <- .calcFactorQuantile(data=x, lib.size=lib.size, p=0.75)
			if( is.null(refColumn) ) refColumn <- which.min(abs(f75-mean(f75)))
			if(length(refColumn)==0 | refColumn < 1 | refColumn > ncol(x)) refColumn <- 1
			f <- rep(NA,ncol(x))
			for(i in 1:ncol(x))
				f[i] <- .calcFactorTMM(obs=x[,i],ref=x[,refColumn], libsize.obs=lib.size[i], libsize.ref=lib.size[refColumn], logratioTrim=logratioTrim, sumTrim=sumTrim, doWeighting=doWeighting, Acutoff=Acutoff)
			f
		},
		TMMwsp = {
			f75 <- .calcFactorQuantile(data=x, lib.size=lib.size, p=0.75)
			if( is.null(refColumn) ) refColumn <- which.min(abs(f75-mean(f75)))
			if(length(refColumn)==0 | refColumn < 1 | refColumn > ncol(x)) refColumn <- 1
			f <- rep(NA,ncol(x))
			for(i in 1:ncol(x))
				f[i] <- .calcFactorTMMwsp(obs=x[,i],ref=x[,refColumn], libsize.obs=lib.size[i], libsize.ref=lib.size[refColumn], logratioTrim=logratioTrim, sumTrim=sumTrim, doWeighting=doWeighting, Acutoff=Acutoff)
			f
		},
		RLE = .calcFactorRLE(x)/lib.size,
		upperquartile = .calcFactorQuantile(x,lib.size,p=p),
		none = rep(1,ncol(x))
	)

#	Factors should multiple to one
	f <- f/exp(mean(log(f)))

#	Output
	f
}

.calcFactorRLE <- function (data)
{
	gm <- exp(rowMeans(log(data)))
	apply(data, 2, function(u) median((u/gm)[gm > 0]))
}

.calcFactorQuantile <- function (data, lib.size, p=0.75)
{
#	i <- apply(data<=0,1,all)
#	if(any(i)) data <- data[!i,,drop=FALSE]
	y <- t(t(data)/lib.size)
	f <- apply(y,2,function(x) quantile(x,p=p))
#	f/exp(mean(log(f)))
}

.calcFactorTMM <- function(obs, ref, libsize.obs=NULL, libsize.ref=NULL, logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10)
#	TMM between two libraries
{
	obs <- as.numeric(obs)
	ref <- as.numeric(ref)

	if( is.null(libsize.obs) ) nO <- sum(obs) else nO <- libsize.obs
	if( is.null(libsize.ref) ) nR <- sum(ref) else nR <- libsize.ref

	logR <- log2((obs/nO)/(ref/nR))          # log ratio of expression, accounting for library size
	absE <- (log2(obs/nO) + log2(ref/nR))/2  # absolute expression
	v <- (nO-obs)/nO/obs + (nR-ref)/nR/ref   # estimated asymptotic variance

#	remove infinite values, cutoff based on A
	fin <- is.finite(logR) & is.finite(absE) & (absE > Acutoff)

	logR <- logR[fin]
	absE <- absE[fin]
	v <- v[fin]

	if(max(abs(logR)) < 1e-6) return(1)

#	taken from the original mean() function
	n <- length(logR)
	loL <- floor(n * logratioTrim) + 1
	hiL <- n + 1 - loL
	loS <- floor(n * sumTrim) + 1
	hiS <- n + 1 - loS

#	keep <- (rank(logR) %in% loL:hiL) & (rank(absE) %in% loS:hiS)
#	a fix from leonardo ivan almonacid cardenas, since rank() can return
#	non-integer values when there are a lot of ties
	keep <- (rank(logR)>=loL & rank(logR)<=hiL) & (rank(absE)>=loS & rank(absE)<=hiS)

	if(doWeighting)
		f <- sum(logR[keep]/v[keep], na.rm=TRUE) / sum(1/v[keep], na.rm=TRUE)
	else
		f <- mean(logR[keep], na.rm=TRUE)

#	Results will be missing if the two libraries share no features with positive counts
#	In this case, return unity
	if(is.na(f)) f <- 0
	2^f
}

.calcFactorTMMwsp <- function(obs, ref, libsize.obs=NULL, libsize.ref=NULL, logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10)
#	TMM with pairing of singleton positive counts between the obs and ref libraries
#	Gordon Smyth
#	Created 19 Sep 2018. Last modified 23 April 2019.
{
	obs <- as.numeric(obs)
	ref <- as.numeric(ref)

#	epsilon serves as floating-point zero
	eps <- 1e-14

#	Identify zero counts
	pos.obs <- (obs > eps)
	pos.ref <- (ref > eps)
	npos <- 2L * pos.obs + pos.ref

#	Remove double zeros and NAs
	i <- which(npos==0L | is.na(npos))
	if(length(i)) {
		obs <- obs[-i]
		ref <- ref[-i]
		npos <- npos[-i]
	}

#	Check library sizes	
	if(is.null(libsize.obs)) libsize.obs <- sum(obs)
	if(is.null(libsize.ref)) libsize.ref <- sum(ref)

#	Pair up as many singleton positives as possible
#	The unpaired singleton positives are discarded so that no zeros remain
	zero.obs <- (npos == 1L)
	zero.ref <- (npos == 2L)
	k <- (zero.obs | zero.ref)
	n.eligible.singles <- min( sum(zero.obs), sum(zero.ref))
	if(n.eligible.singles > 0L) {
		refk <- sort(ref[k],decreasing=TRUE)[1:n.eligible.singles]
		obsk <- sort(obs[k],decreasing=TRUE)[1:n.eligible.singles]
		obs <- c(obs[!k],obsk)
		ref <- c(ref[!k],refk)
	} else {
		obs <- obs[!k]
		ref <- ref[!k]
	}

#	Any left?
	n <- length(obs)
	if(n==0L) return(1)

#	Compute M and A values
	obs.p <- obs / libsize.obs
	ref.p <- ref / libsize.ref
	M <- log2( obs.p / ref.p )
	A <- 0.5 * log2( obs.p * ref.p )

#	If M all zero, return 1
	if(max(abs(M)) < 1e-6) return(1)

#	M order, breaking ties by shrunk M
	obs.p.shrunk <- (obs+0.5) / (libsize.obs+0.5)
	ref.p.shrunk <- (ref+0.5) / (libsize.ref+0.5)
	M.shrunk <- log2( obs.p.shrunk / ref.p.shrunk )
	o.M <- order(M, M.shrunk)

#	A order
	o.A <- order(A)

#	Trim
	loM <- as.integer(n * logratioTrim) + 1L
	hiM <- n + 1L - loM
	keep.M <- rep(FALSE,n)
	keep.M[o.M[loM:hiM]] <- TRUE
	loA <- as.integer(n * sumTrim) + 1L
	hiA <- n + 1L - loA
	keep.A <- rep(FALSE,n)
	keep.A[o.A[loA:hiA]] <- TRUE
	keep <- keep.M & keep.A
	M <- M[keep]

#	Average the M values
	if(doWeighting) {
		obs.p <- obs.p[keep]
		ref.p <- ref.p[keep]
		v <- (1-obs.p)/obs.p/libsize.obs + (1-ref.p)/ref.p/libsize.ref
		w <- (1+1e-6) / (v+1e-6)
		TMM <- sum(w*M) / sum(w)
	} else {
		TMM <- mean(M)
	}

	2^TMM
}
