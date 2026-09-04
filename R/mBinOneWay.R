mbinOneWay <- function(y, z, design=NULL, group=NULL, offset=NULL, weights=NULL, maxit=50, tol=1e-8, nthreads=1L)
#	Fit one-way layout binomial models using weighted mean
#	Created by Lizhong Chen
#	11 Apr 2025.  Last modified 18 Apr 2025.
{
#	check counts
	y <- as.matrix(y)
	z <- as.matrix(z)

	ngenes <- nrow(y)
	nlibs  <- ncol(y)

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

#	If necessary, the group factor is computed from the design matrix.
#	However, if group is supplied, we can avoid creating a design matrix altogether.
	if(is.null(group)) {	
		if(is.null(design)) {
			group <- factor(rep_len(1L,nlibs))
		} else {
			design <- as.matrix(design)
			if(nrow(design) != nlibs) stop("nrow(design) disagrees with ncol(y)")
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
		designunique  <- design[firstjofgroup,,drop=FALSE]
		
#		It is just a group indicator matrix?
		if(sum(designunique==1)==ngroups && sum(designunique==0)==(ngroups-1L)*ngroups) design <- NULL
	}

#	Cycle through groups
	mu   <- prop
	beta <- dev <- matrix(0,ngenes,ngroups)

	for (g in seq_len(ngroups)) {
		j   <- which(i==g)
		out <- .Call(.cxx_bin_one_group, prop[,j,drop=FALSE], offset[,j,drop=FALSE], weights[,j,drop=FALSE], maxit, tol, nthreads)
		beta[,g] <- out$coef
		dev[,g]  <- out$deviance
		mu[,j]   <- out$fitted.values
	}
	deviance <- rowSums(dev)

#	If necessary, reformat the beta's to reflect the original design.
	if(!is.null(design)) {
		beta <- t(solve(designunique,t(beta)))
		rownames(beta) <- rownames(y)
	} else {
		dimnames(beta) <- list(rownames(y),levg)
	}

	list(coefficients=beta, fitted.values=mu, deviance=deviance)
}
