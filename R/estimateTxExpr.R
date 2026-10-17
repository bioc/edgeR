estimateTxExpr <- function(lengths,counts,compatibility,maxit=1000L,eps=1e-6,start.rpb=NULL)
#	Estimate transcript expressions from equivalence class counts and compatibility matrix.
#	Gordon Smyth
#	Created 9 Oct 2022. Last revised 16 Oct 2026.
{
#	Number of transcripts
	ntx <- length(lengths)

#	Coerce counts to matrix, allowing for multiple input vectors from bootstraps or simulations
	counts <- as.matrix(counts)

#	Number of equivalence classes
	nec <- nrow(counts)

#	Number of input profiles
	nsim <- ncol(counts)

#	Starting values
	if(is.null(start.rpb)) {
		A <- matrix(1,ntx,nsim)
	} else {
		A <- matrix(start.rpb,ntx,nsim)
	}
	row.names(A) <- names(lengths)
	Y <- A

#	Form compatibility matrix
	if(is.matrix(compatibility)) {
		X <- compatibility
	} else {
		X <- matrix(0,nec,ntx)
		Row <- rep(1:nec,lengths(compatibility))
		Col <- unlist(compatibility)
		X[Row + (Col-1)*nec] <- 1
	}

#	EM algorithm
	iter <- 0L
	maxdiff <- 1
	Zero <- (counts==0)
	repeat {
		iter <- iter + 1L
		Aold <- A
		maxdiffold <- maxdiff
		R <- counts / (X %*% A)
		R[Zero] <- 0
		Y <- A * crossprod(X,R)
		A <- Y / lengths
		maxdiff <- max(abs(A-Aold))
		contraction.factor <- maxdiff/maxdiffold
		delta <- maxdiff/max(A) * contraction.factor/(1-contraction.factor)
		if(delta < eps) break
		if(iter > maxit) break
	}

#	Dimnames
	TxNames <- names(lengths)
	if(is.null(TxNames)) TxNames <- paste0("Tx",1:ntx)
	rownames(A) <- rownames(Y) <- TxNames
	if(nsim > 1L) {
		SimNames <- colnames(counts)
		if(is.null(SimNames)) SimNames <- paste0("S",1:nsim)
		colnames(A) <- colnames(Y) <- SimNames
	}

	list(rpb=drop(A),est.count=drop(Y),niter=iter,delta=delta,contraction.factor=contraction.factor)
}
