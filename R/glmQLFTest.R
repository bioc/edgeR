#  FIT QUASI-LIKELIHOOD GENERALIZED LINEAR MODELS

glmQLFit <- function(y, ...)
UseMethod("glmQLFit")

glmQLFit.DGEList <- function(y, design=NULL, dispersion=NULL, abundance.trend=TRUE, robust=FALSE, winsor.tail.p=c(0.05, 0.1), legacy=FALSE, top.proportion=0.1, keep.unit.mat=FALSE,...)
# 	Fit NB GLMs and estimate QL dispersions with empirical Bayes moderation.
# 	Yunshun Chen, Aaron Lun, Lizhong Chen, Gordon Smyth
#	Created 5 November 2014. Last modified 8 Apr 2024.
{
#	The design matrix defaults to the oneway layout defined by y$samples$group.
#	If there is only one group, then the design matrix is left NULL so that a
#	matrix with a single intercept column will be set later by glmFit.default.
	if(is.null(design)) {
		design <- y$design
		if(is.null(design)) {
			group <- droplevels(as.factor(y$samples$group))
			if(nlevels(group) > 1L) design <- model.matrix(~y$samples$group)
		}
	}

	if(is.null(y$AveLogCPM)) y$AveLogCPM <- aveLogCPM(y)

	if(legacy && is.null(dispersion)) {
		dispersion <- y$trended.dispersion
		if(is.null(dispersion)) dispersion <- y$common.dispersion
		if(is.null(dispersion)) stop("No dispersion values found in DGEList object.")
	}	

	offset <- getOffset(y)

	fit <- glmQLFit.default(y=y$counts, design=design, dispersion=dispersion, offset=offset, lib.size=NULL, abundance.trend=abundance.trend, AveLogCPM=y$AveLogCPM, robust=robust, winsor.tail.p=winsor.tail.p, weights=y$weights, legacy=legacy, top.proportion=top.proportion, keep.unit.mat=keep.unit.mat, ...)

	fit$samples <- y$samples
	fit$genes <- y$genes
	fit$AveLogCPM <- y$AveLogCPM
	fit
}

glmQLFit.SummarizedExperiment <- function(y, design=NULL, dispersion=NULL, abundance.trend=TRUE, robust=FALSE, winsor.tail.p=c(0.05, 0.1), legacy=FALSE, top.proportion=0.1,keep.unit.mat=FALSE,...)
#	Created 3 April 2020. Last modified 8 April 2023.
{
	y <- SE2DGEList(y)
	glmQLFit.DGEList(y, design=design, dispersion=dispersion, abundance.trend=abundance.trend, robust=robust, winsor.tail.p=winsor.tail.p, legacy=legacy, top.proportion=top.proportion,keep.unit.mat=keep.unit.mat,...)
}

glmQLFit.default <- function(y, design=NULL, dispersion=NULL, offset=NULL, lib.size=NULL, weights=NULL, 
	abundance.trend=TRUE, AveLogCPM=NULL, covariate.trend=NULL,
	robust=FALSE, winsor.tail.p=c(0.05, 0.1),
	legacy=FALSE, top.proportion=0.1, keep.unit.mat=FALSE, ...)
# 	Fits genewise GLMs and estimates quasi-likelihood dispersions with empirical Bayes moderation.
# 	Originally part of glmQLFTest created by Davis McCarthy and Gordon Smyth, 13 Jan 2012.
#	DF adjustment for zeros added by Aaron Lun and Gordon Smyth, 7 Jan 2014.
#	Split from glmQLFTest as a separate function by Aaron Lun and Yunshun Chen, 15 Sep 2014.
#	Bias adjustment for deviance and DF added by Lizhong Chen and Gordon Smyth, 8 Nov 2022.
#	C++ replaced with pure C by Lizhong Chen, 6 May 2024.
#	legacy=FALSE argument passed to squeezeVar(), 1 Aug 2024.
#	Last modified 8 Oct 2024.
{
#	Check y
	y <- as.matrix(y)

#	Check design
	if(is.null(design)) {
		design <- matrix(1,ncol(y),1)
		rownames(design) <- colnames(y)
		colnames(design) <- "Intercept"
	}

#	Check AveLogCPM
	if(is.null(AveLogCPM)) AveLogCPM <- aveLogCPM(y, offset=offset, lib.size=lib.size, weights=weights, dispersion=dispersion)

#	Check weights
	weights <- .compressWeights(y,weights)

#	Check offsets
	offset <- .compressOffsets(y,offset,lib.size)

#	Check dispersion
	if(is.null(dispersion)) {
		if(legacy){
			stop("No dispersion values provided.")
		} else {
			if(top.proportion < 0 || top.proportion > 1) stop("top.proportion should be between 0 and 1.")
			i <- (AveLogCPM >= quantile(AveLogCPM, probs=1-top.proportion))
			dispersion <- estimateGLMCommonDisp(y[i,,drop=FALSE], design=design, offset=offset[i,,drop=FALSE], weights=weights[i,,drop=FALSE])
		}
	}	

	fit <- glmFit.default(y, design=design, dispersion=dispersion, offset=offset, lib.size=lib.size, weights=weights, ...)

#	Store AveLogCPM for the computation of average ql dispersion
	AveLogCPM2 <- AveLogCPM

#	Covariate for trended prior for quasi-dispersion
	if(is.null(covariate.trend)) {
		if(abundance.trend) {
			fit$AveLogCPM <- AveLogCPM
		} else {
			AveLogCPM <- NULL
		}
	} else{
		AveLogCPM <- covariate.trend
	}

#	Setting the residual deviances and df
	if(legacy) {

#		Old-style adjustment retained as option for backward compatibility.
#		Adjust df.residual for fitted values at zero.
		zerofit <- (fit$fitted.values < 1e-4) & (fit$counts < 1e-4)
		df.residual <- .residDF(zerofit, fit$design)
		fit$df.residual.zeros <- df.residual
		s2 <- fit$deviance / df.residual
		s2[df.residual==0L] <- 0

	} else {

#		New-style adjustment.
#		Deviance and df.residual are both adjusted for bias and variance.

#		Check dispersion for CompressedMatrix
		dispersion.mat <- .compressDispersions(y,dispersion)

#		Compute average quasi dispersion
		ave.ql.disp  <- .Call(.cxx_compute_ave_qd,y,fit$fitted.values,design,dispersion.mat,AveLogCPM2,weights)

# 		Refit using the scaled dispersion by average quasi dispersion 
		fit <- glmFit(y, design=design, dispersion=dispersion/ave.ql.disp, offset=offset, lib.size=lib.size, weights=weights, ...)
		fit$dispersion <- dispersion

#		Prepare outputs

#		First the adjusted deviance, df and estimated quasi-dispersions
		out <- .Call(.cxx_compute_adj_vec,y,fit$fitted.values,design,dispersion.mat,ave.ql.disp,weights)
		s2  <- out$s2
		fit$df.residual.adj <- df.residual <- out$df
		fit$deviance.adj	<- out$deviance
		fit$average.ql.dispersion		  <- ave.ql.disp

#		Second the adjusted unit deviance, df and leverage matrix
		if(keep.unit.mat){
			out <- .Call(.cxx_compute_adj_mat,y,fit$fitted.values,design,dispersion.mat,ave.ql.disp,weights)
			fit$leverage		  <- out$leverage
			fit$unit.deviance.adj <- out$unit.deviance
			fit$unit.df.adj	   <- out$unit.df
		}
	}

#	Empirical Bayes moderation of quasi-likelihood dispersions
	s2.fit <- squeezeVar(s2,df=df.residual,covariate=AveLogCPM,robust=robust,winsor.tail.p=winsor.tail.p,legacy=legacy)

#	Storing results
	fit$df.prior <- s2.fit$df.prior
	fit$s2.post  <- s2.fit$var.post
	fit$s2.prior <- s2.fit$var.prior
	fit
}

glmQLFTest <- function(glmfit, coef=ncol(glmfit$design), contrast=NULL, poisson.bound=TRUE)
#	Quasi-likelihood F-tests for DGE glms.
#	Davis McCarthy, Gordon Smyth, Aaron Lun.
#	Created 18 Feb 2011. Last modified 1 Aug 2024.
{
	if(!is(glmfit,"DGEGLM")) stop("glmfit must be an DGEGLM object produced by glmQLFit") 
	if(is.null(glmfit$s2.post)) stop("need to run glmQLFit before glmQLFTest")

#	add check in glmLRT function for new QL method
#	fitting null model using working dispersion
	out <- glmLRT(glmfit, coef=coef, contrast=contrast)

#	Older code stored df adjusted for exact zeros in df.residual.zero.
#	New code stores adjusted df in df.residual.adj.
	if(is.null(glmfit$df.residual.zeros)) {
		df.residual <- glmfit$df.residual.adj
		poisson.bound <- FALSE
	} else {
		df.residual <- glmfit$df.residual.zeros
	}

#	Compute the QL F-statistic
	F.stat <- out$table$LR / out$df.test / glmfit$s2.post
	df.total <- glmfit$df.prior + df.residual
	df.residual.total <- sum(glmfit$df.residual)
	df.total <- pmin(df.total, df.residual.total)

#	Compute p-values from the QL F-statistic
	F.pvalue <- pf(F.stat, df1=out$df.test, df2=df.total, lower.tail=FALSE, log.p=FALSE)

#	Ensure is not more significant than chisquare test with Poisson variance.
#	This step is obsolete in the new code.
	if(poisson.bound) {
		i <- .isBelowPoissonBound(glmfit)
		if(any(i)) {
			pois.fit <- glmfit[i,]
			pois.fit <- glmFit(pois.fit$counts, design=pois.fit$design, offset=pois.fit$offset, weights=pois.fit$weights, start=pois.fit$unshrunk.coefficients, dispersion=0)
			pois.res <- glmLRT(pois.fit, coef=coef, contrast=contrast) 
			F.pvalue[i] <- pmax(F.pvalue[i], pois.res$table$PValue)
		}
	}

	out$table$LR <- out$table$PValue <- NULL
	out$table$F <- F.stat
	out$table$PValue <- F.pvalue
	out$df.total <- df.total

	out
}

.isBelowPoissonBound <- function(glmfit) 
# A convenience function to avoid generating temporary matrices.
{
	disp <- makeCompressedMatrix(glmfit$dispersion, dim(glmfit$counts), byrow=FALSE)
	s2 <- makeCompressedMatrix(glmfit$s2.post, dim(glmfit$counts), byrow=FALSE)
	out <- .Call(.cxx_check_poisson_bound, glmfit$fitted.values, disp, s2)
	return(out)
}

plotQLDisp <- function(glmfit, xlab="Average Log2 CPM", ylab="Quarter-Root Mean Deviance", pch=16, cex=0.2, col.shrunk="red", col.trend="blue", col.raw="black", ...)
# 	Plots raw and empirical Bayes moderated quasi dispersion estimates.
# 	Originally part of glmQLFTest created by Davis McCarthy and Gordon Smyth, 13 Jan 2012.
#	DF adjustment for zeros added by Aaron Lun and Gordon Smyth, 7 Jan 2014.
#	Split from glmQLFTest as separate function by Aaron Lun and Yunshun Chen, 15 Sep 2014.
#	Bias adjustment for deviance and DF added by Lizhong Chen and Gordon Smyth, 8 Nov 2022.
#	Last modified 22 Jan 2023.
{
	if(is.null(glmfit$s2.post)) stop("need to run glmQLFit before plotQLDisp")

#	Make sure average logCPM is available
	A <- glmfit$AveLogCPM
	if(is.null(A)) A <- aveLogCPM(glmfit)

#	Older code put df adjusted for exact zeros in df.residual.zero.
#	Newer code puts adjusted df in df.residual.adj.
	if(is.null(glmfit$df.residual.zeros)) {
		df.residual <- glmfit$df.residual.adj
		deviance <- glmfit$deviance.adj
	} else {
		df.residual <- glmfit$df.residual.zeros
		deviance <- glmfit$deviance
	}
	s2 <- deviance / df.residual
	s2[df.residual < 1e-8] <- 0

	plot(A, sqrt(sqrt(s2)),xlab=xlab, ylab=ylab, pch=pch, cex=cex, col=col.raw, ...)
	points(A, sqrt(sqrt(glmfit$s2.post)), pch=pch, cex=cex, col=col.shrunk)
	if(identical(length(glmfit$s2.prior),1L)) { 
		abline(h=sqrt(sqrt(glmfit$s2.prior)), col=col.trend)
	} else {
		o <- order(A)
		lines(A[o], sqrt(sqrt(glmfit$s2.prior[o])), col=col.trend, lwd=2)
	}
	legend("topright", lty=c(-1,-1,1), pch=c(pch,pch,-1), col=c(col.raw,col.shrunk,col.trend), pt.cex=0.7, lwd=2, legend=c("Raw","Squeezed", "Trend"))

	invisible(list(x=A,y=sqrt(sqrt(s2))))
}
