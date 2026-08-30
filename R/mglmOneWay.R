designAsFactor <- function(design)
#	Construct a factor from the unique rows of a matrix
#	Gordon Smyth
#	11 March 2011.  Last modified 19 March 2011.
{
	design <- as.matrix(design)
	z <- (exp(1)+pi)/5
	g <- factor(rowMeans(design*z^(col(design)-1)))
	levels(g) <- seq_len(length(levels(g)))
	g
}

mglmOneWay <- function(y,design=NULL,group=NULL,dispersion=0,offset=0,weights=NULL,coef.start=NULL,maxit=50,tol=1e-10,nthreads=1L)
#	Fit multiple negative binomial glms with log link
#	by Fisher scoring with
#	only a single explanatory factor in the model
#	Gordon Smyth
#	Aapted to use .compress and C functions by Aaron Lun
#	11 March 2011.  Last modified 24 August 2017.
{
	y      <- as.matrix(y)
	ngenes <- nrow(y)
	nlibs  <- ncol(y)

	offset     <- .compressOffsets(y, offset=offset)
	dispersion <- .compressDispersions(y, dispersion)
	weights    <- .compressWeights(y, weights)

#	If necessary, the group factor is computed from the design matrix.
#	However, if group is supplied, we can avoid creating a design matrix altogether.
	if(is.null(group)) {	
		if(is.null(design)) {
			group <- factor(rep_len(1L,nlibs))
		} else {
			design <- as.matrix(design)
			group  <- designAsFactor(design)
		}
	} else {
		group <- as.factor(group)
	}

#	Convert factor to integer levels for efficiency
	levg    <- levels(group)
	ngroups <- length(levg)
	i       <- as.integer(group)

	if(!is.null(design)) {
		if(ncol(design)!=ngroups) stop("design matrix is not equivalent to a oneway layout")

#		Reduce to representative design matrix, based on the column in which each group appears first.
		firstjofgroup <- match(levg, group)
		designunique <- design[firstjofgroup,,drop=FALSE]

#		It is just a group indicator matrix?
		if(sum(designunique==1)==ngroups && sum(designunique==0)==(ngroups-1L)*ngroups) design <- NULL

#		If necessary, convert starting values to group fitted values
		if(!is.null(design) && !is.null(coef.start)) coef.start <- coef.start %*% t(designunique)
	}

	if (is.null(coef.start)) {
		coef.start.mat <- matrix(NA_real_, ngenes, ngroups)
	} else {
		coef.start.mat <- as.matrix(coef.start)
		if (!is.double(coef.start.mat)) storage.mode(coef.start.mat) <- "double"
	}

#	Cycle through groups on C level
	output   <- .Call(.cxx_fit_one_way, y, offset, dispersion, weights, i - 1L, ngroups, maxit, tol, coef.start.mat, nthreads)
	beta     <- output$coefficients

#	Reset -Inf values to finite value to simplify calculations downstream
#	Now redundant: the C fitter (glm_one_group_vec / fit_one_way_mat) already returns -1e8 for
#	all-zero groups, so this pmax is a no-op. Kept commented for reference.
#	beta <- pmax(beta,-1e8)

#	Fitted values from group-wise beta's
	mu <- .Call(.cxx_get_one_way_fitted, beta, offset, i-1L, nthreads)
	dimnames(mu) <- dimnames(y)

#	If necessary, reformat the beta's to reflect the original design.
	if(!is.null(design)) {
		beta <- t(solve(designunique,t(beta)))
		rownames(beta) <- rownames(y)
	} else {
		dimnames(beta) <- list(rownames(y),levg)
	}

	list(coefficients=beta,fitted.values=mu)
}
