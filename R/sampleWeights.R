sampleWeights <- function(glmfit)
#	Empirical sample weights from adjusted unit deviances
#	assuming that the quasi-dispersion of each observation is s2_g / w_i.
#	Handles both DGEGLM (NB, from glmQLFit) and DGEBIN (binomial, from binQLFit)
#	fits, with nthreads read from the fit object.
#	Lizhong Chen, Gordon Smyth
#	Created 5 Aug 2024. Last modified 16 Jul 2026.
{
#	shared setup
	weights  <- .compressWeights(glmfit$counts,glmfit$weights)
	design   <- glmfit$design
	nthreads <- glmfit$nthreads

	if(is(glmfit,"DGEGLM")) {
#		check new QL pipeline
		if(is.null(glmfit$average.ql.dispersion)) stop("Please run glmQLFit with legacy=FALSE first")

#		negative-binomial pipeline. The adjusted unit deviance/df and the
#		relative sample-weight reduction (using the prior s2, floored at 1) are
#		folded into a single C call returning the per-sample weights.
		dispersion.mat <- .compressDispersions(glmfit$counts,glmfit$dispersion)

#		out <- .Call(.cxx_compute_adj_mat,glmfit$counts,glmfit$fitted.values,design,dispersion.mat,glmfit$average.ql.dispersion,weights,nthreads)
#		s2 <- pmax(glmfit$s2.prior,1)
#		w  <- colSums(out$unit.deviance * as.matrix(weights) / s2) / colSums(out$unit.df)
#		w  <- log(w)
#		w  <- exp(mean(w)-w)
		w <- .Call(.cxx_sample_weights,glmfit$counts,glmfit$fitted.values,design,dispersion.mat,glmfit$average.ql.dispersion,weights,glmfit$s2.prior,nthreads)

	} else if(is(glmfit,"DGEBIN")) {
#		binomial pipeline. The C call takes the raw counts and counts2; coverage,
#		proportion and scaled weights are derived inside C, and the sample-weight
#		reduction (s2=1) is folded into the same call.

#		coverage <- glmfit$counts + glmfit$counts2
#		weights0 <- as.matrix(weights) * coverage
#		prop     <- glmfit$counts / pmax(coverage,1)
#		out <- .Call(.cxx_compute_adj_mat_bin,prop,glmfit$fitted.values,design,coverage,weights0,nthreads)
#		w  <- colSums(out$unit.deviance * as.matrix(weights)) / colSums(out$unit.df)
#		w  <- log(w)
#		w  <- exp(mean(w)-w)
		w <- .Call(.cxx_sample_weights_bin,glmfit$counts,glmfit$counts2,glmfit$fitted.values,design,weights,nthreads)

	} else {
		stop("glmfit must be a DGEGLM or DGEBIN object")
	}

#	return the per-sample weights as a compressedMatrix repeated across genes
	makeCompressedMatrix(w, dim(glmfit$counts), byrow=TRUE)
}
