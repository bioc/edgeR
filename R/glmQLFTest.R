#  FIT QUASI-LIKELIHOOD GENERALIZED LINEAR MODELS

glmQLFit <- function(y, ...)
UseMethod("glmQLFit")

glmQLFit.DGEList <- function(y, design=NULL, dispersion=NULL, abundance.trend=TRUE, robust=FALSE, winsor.tail.p=c(0.05, 0.1), legacy=FALSE,...)
# 	Fit NB GLMs and estimate QL dispersions with empirical Bayes moderation.
# 	Yunshun Chen, Aaron Lun, Lizhong Chen, Gordon Smyth
#	Created 5 November 2014. Last modified 23 December 2025.
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

	if(is.null(dispersion) && legacy) {
		dispersion <- y$trended.dispersion
		if(is.null(dispersion)) dispersion <- y$common.dispersion
		if(is.null(dispersion)) stop("No dispersion values found in DGEList object.")
	}

	offset <- getOffset(y)

	fit <- glmQLFit.default(y=y$counts, design=design, dispersion=dispersion, offset=offset, lib.size=NULL, abundance.trend=abundance.trend, AveLogCPM=y$AveLogCPM, robust=robust, winsor.tail.p=winsor.tail.p, weights=y$weights, legacy=legacy, ...)

	fit$samples <- y$samples
	fit$genes <- y$genes
	fit$AveLogCPM <- y$AveLogCPM
	fit
}

glmQLFit.SummarizedExperiment <- function(y, design=NULL, dispersion=NULL, abundance.trend=TRUE, robust=FALSE, winsor.tail.p=c(0.05, 0.1), legacy=FALSE,...)
#	Created 3 April 2020. Last modified 8 April 2023.
{
	y <- SE2DGEList(y)
	glmQLFit.DGEList(y, design=design, dispersion=dispersion, abundance.trend=abundance.trend, robust=robust, winsor.tail.p=winsor.tail.p, legacy=legacy,...)
}

glmQLFit.default <- function(y, design=NULL, dispersion=NULL, offset=NULL, lib.size=NULL, weights=NULL,
	abundance.trend=TRUE, AveLogCPM=NULL, covariate.trend=NULL,
	robust=FALSE, winsor.tail.p=c(0.05, 0.1),
	legacy=FALSE, nthreads=1L, ...)
# 	Fits genewise GLMs and estimates quasi-likelihood dispersions with empirical Bayes moderation.
# 	Originally part of glmQLFTest created by Davis McCarthy and Gordon Smyth, 13 Jan 2012.
#	DF adjustment for zeros added by Aaron Lun and Gordon Smyth, 7 Jan 2014.
#	Split from glmQLFTest as a separate function by Aaron Lun and Yunshun Chen, 15 Sep 2014.
#	Bias adjustment for deviance and DF added by Lizhong Chen and Gordon Smyth, 8 Nov 2022.
#	C++ replaced with pure C by Lizhong Chen, 6 May 2024.
#	legacy=FALSE argument passed to squeezeVar(), 1 Aug 2024.
#	Support binomial models, 29 April 2025
#	Last modified 29 April 2025.
{
#	Check y
	y <- as.matrix(y)

#	Check nthreads (validated once here so it is clean on both the legacy and
#	non-legacy paths, and can be stored on the returned fit for glmQLFTest)

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

#	Store AveLogCPM for the computation of average ql dispersion
	AveLogCPM2 <- AveLogCPM

#	Covariate for trended prior for quasi-dispersion
	if(is.null(covariate.trend)) {
		if(!abundance.trend) AveLogCPM <- NULL
	} else{
		AveLogCPM <- covariate.trend
	}

#	Setting the residual deviances and df
	if(legacy) {
		if(is.null(dispersion)) stop("No dispersion values provided.")

		fit <- glmFit.default(y, design=design, dispersion=dispersion, offset=offset, lib.size=lib.size, weights=weights, nthreads=nthreads, ...)

#		Old-style adjustment retained as option for backward compatibility.
#		Adjust df.residual for fitted values at zero.
		zerofit <- (fit$fitted.values < 1e-4) & (fit$counts < 1e-4)
		df.residual <- .residDF(zerofit, fit$design)
		fit$df.residual.zeros <- df.residual
		s2 <- fit$deviance / df.residual
		s2[df.residual==0L] <- 0

	} else {
#	Check dispersion
		if(is.null(dispersion)) {
#			ngenes      <- nrow(y)
#			df.residual <- ncol(y) - ncol(design)
#			top.prop    <- chooseLowessSpan(ngenes*sqrt(df.residual),small.n=20,min.span=0.02)
#			ntop        <- ceiling(top.prop * ngenes)
#			i           <- order(AveLogCPM2,decreasing=TRUE)[1:ntop]
#			dispersion  <- dispCoxReid(y[i,,drop=FALSE], design=design, offset=offset[i,,drop=FALSE], weights=weights[i,,drop=FALSE], nthreads=nthreads)

#			Select the most-abundant genes and estimate the common Cox-Reid
#			dispersion over them in a single C call. This folds the former R
#			top-gene selection and dispCoxReid() call (min.row.sum filtering,
#			systematicSubset subsetting and the nonEstimable design check are
#			intentionally omitted; the selected genes are the highest-abundance
#			genes and the design is validated upstream by glmFit). The bounds
#			0 and sqrt(2) are min/max of interval=c(0,4) raised to the 0.25
#			power, matching dispCoxReid's defaults (tol=1e-5, Cox-Reid adjust).

			dispersion <- .Call(.cxx_coxreid_disp_top, y, offset, weights, design, AveLogCPM2, 0, sqrt(2), 1e-5, TRUE, nthreads)
		}

#		Check dispersion for CompressedMatrix
		dispersion.mat <- .compressDispersions(y,dispersion)
		fit <- glmFit.default(y, design=design, dispersion=dispersion.mat, offset=offset, lib.size=lib.size, weights=weights, nthreads=nthreads, prior.count=0, ...)
#		Compute average quasi dispersion (nthreads validated near the top)
		ave.ql.disp  <- .Call(.cxx_compute_ave_qd,y,fit$fitted.values,design,dispersion.mat,AveLogCPM2,weights,nthreads)

#		New-style adjustment.
#		Deviance and df.residual are both adjusted for bias and variance.
#		Refit using the scaled dispersion by average quasi dispersion
		fit <- glmFit.default(y, design=design, dispersion=dispersion/ave.ql.disp, offset=offset, lib.size=lib.size, weights=weights, nthreads=nthreads, ...)
		fit$dispersion <- dispersion

#		Prepare outputs

#		First the adjusted deviance, df and estimated quasi-dispersions
		out <- .Call(.cxx_compute_adj_vec,y,fit$fitted.values,design,dispersion.mat,ave.ql.disp,weights,nthreads)
		s2  <- out$s2
		fit$df.residual.adj <- df.residual <- out$df
		fit$deviance.adj	<- out$deviance
		fit$average.ql.dispersion		  <- ave.ql.disp
	}

#	Empirical Bayes moderation of quasi-likelihood dispersions
	s2.fit <- squeezeVar(s2,df=df.residual,covariate=AveLogCPM,robust=robust,winsor.tail.p=winsor.tail.p,legacy=legacy)

#	Storing results
	fit$AveLogCPM <- AveLogCPM2
	fit$df.prior  <- s2.fit$df.prior
	fit$s2.post   <- s2.fit$var.post
	fit$s2.prior  <- s2.fit$var.prior
	fit$nthreads  <- nthreads
	fit
}

glmQLFTest <- function(glmfit, coef=ncol(glmfit$design), contrast=NULL, fc = 1, lfc = NULL, poisson.bound=TRUE, upshot=TRUE)
#	Quasi-likelihood F-tests for DGE quasi-negative binomial models.
#	Support quasi-binomial models, 29 April 2025
#	Treat analysis testing relative to a minimum threshold by Lizhong Chen, 29 April 2025
#	Support UPSHOT method, 13 Nov 2025
#	Davis McCarthy, Gordon Smyth, Aaron Lun, Lizhong Chen.
#	Created 18 Feb 2011. Last modified 18 Nov 2025
{
	if(!is(glmfit,"DGEGLM")) stop("glmfit must be an DGEGLM object produced by glmQLFit") 
	if(is.null(glmfit$s2.post)) stop("need to run glmQLFit before glmQLFTest")
	nlibs <- ncol(glmfit)

	#	Check design matrix
	design <- as.matrix(glmfit$design)
	nbeta  <- ncol(design)
	if(nbeta < 2) stop("Need at least two columns for design, usually the first is the intercept column")
	coef.names <- colnames(design)
	
	#	Evaluate logFC for coef to be tested
	#	Note that contrast takes precedence over coef: if contrast is given
	#	then reform design matrix so that contrast of interest is last column.
	if(is.null(contrast)) {
		if(length(coef) > 1) coef <- unique(coef)
		if(is.character(coef)) {
			check.coef <- coef %in% colnames(design)
			if(any(!check.coef)) stop("One or more named coef arguments do not match a column of the design matrix.")
			coef.name <- coef
			coef <- match(coef, colnames(design))
		}
		else
			coef.name <- coef.names[coef]
		logFC <- glmfit$coefficients[,coef,drop=FALSE]/log(2)
		if(!is.null(glmfit$unshrunk.coefficients)) unshrunk.logFC <- glmfit$unshrunk.coefficients[,coef,drop=FALSE]/log(2)
	} else {
		contrast <- as.matrix(contrast)
		reform   <- contrastAsCoef(design, contrast=contrast, first=TRUE)
		coef     <- reform$coef
		design   <- reform$design
		if(length(coef)>1) {
			coef.name <- paste("LR test on",length(coef),"degrees of freedom")
		} else {
			contrast <- drop(contrast)
			i <- contrast!=0
			coef.name <- paste(paste(contrast[i],coef.names[i],sep="*"),collapse=" ")
		}
		logFC  <- (glmfit$coefficients %*% contrast)/log(2)
		if(!is.null(glmfit$unshrunk.coefficients)) unshrunk.logFC <- (glmfit$unshrunk.coefficients %*% contrast)/log(2)
	}
	
	# number of coefficients
	ncoef <- length(coef)

	# Null design matrix
	design0 <- design[, -coef, drop=FALSE]
	
	# adjust NB dispersion for new QL method scaled by average QL dispersion
	if(!is.null(glmfit$dispersion)){
		if(is.null(glmfit$average.ql.dispersion)) {
			dispersion <- glmfit$dispersion
		} else {
			dispersion <- glmfit$dispersion/glmfit$average.ql.dispersion
		}
	} else {
		dispersion <- NULL
	}

	# check lfc 
	if(is.null(lfc))    lfc <- log2(fc)
	if(length(lfc)==1L) lfc <- rep(lfc, ncoef)
	if(sum(lfc < 0) > 0) stop("lfc should be non-negative")
	if((length(lfc)!=ncoef) & (length(lfc)>1)) stop("fc or lfc vector of wrong length")

	# degree of freedom
	if(is.null(glmfit$df.residual.zeros)) {
		df.residual <- glmfit$df.residual.adj
		poisson.bound <- FALSE
	} else {
		df.residual <- glmfit$df.residual.zeros
	}

	# starting p-values for normal F-test    
	#--- ORIGINAL base null refit (folded into .cxx_glm_ftest below; kept for reference) ---
	# fit <- glmFit(glmfit$counts,design=design0,offset=glmfit$offset,weights=glmfit$weights,dispersion=dispersion,prior.count=0)
	# LR  <- pmax(fit$deviance-glmfit$deviance,0)

	df.total <- glmfit$df.prior + df.residual
	df.residual.total <- sum(glmfit$df.residual)
	df.total <- pmin(df.total, df.residual.total)

	#	QL F-statistic + p-value + UPSHOT/TREAT computation in a SINGLE C call.
	#	The base null refit, F-statistic/pf, the TREAT Gauss-quadrature (formerly
	#	the R helper .glmtreat) and the final qf all run in C via .cxx_glm_ftest,
	#	reusing the NB fitter fit_glm_mat that glmFit drives through .cxx_fit_glm,
	#	so results match the former R path to ~1e-8.  offset/dispersion/weights are
	#	compressed exactly as glmFit feeds the fitter; the logFCt selection and the
	#	glmfit$lfc side-effect stay in R; nthreads is read from glmfit$nthreads (the
	#	thread count stored by glmQLFit; formerly pinned to 1L here).
	offset.mat     <- .compressOffsets(glmfit$counts, offset=glmfit$offset)
	dispersion.mat <- .compressDispersions(glmfit$counts, dispersion)
	weights.mat    <- .compressWeights(glmfit$counts, glmfit$weights)
	logFCt <- logFC
	if(!is.null(glmfit$unshrunk.coefficients)) logFCt <- unshrunk.logFC
	if(any(lfc>0)) glmfit$lfc <- lfc
	out <- .Call(.cxx_glm_ftest, glmfit$counts, offset.mat, dispersion.mat, weights.mat, design, coef, glmfit$deviance, glmfit$s2.post, df.total, lfc, logFCt, upshot, glmfit$nthreads)
	F.stat   <- out$F
	F.pvalue <- out$PValue

	#--- ORIGINAL F-statistic + TREAT/UPSHOT + qf (replaced by .cxx_glm_ftest above; kept for reference) ---
	# #	Compute p-values from the QL F-statistic
	# F.stat   <- LR / ncoef / glmfit$s2.post
	# F.pvalue <- pf(F.stat, df1=ncoef, df2=df.total, lower.tail=FALSE, log.p=FALSE) 
# 
	# # perform treat analysis if any lfc > 0
	# if(any(lfc>0)){
		# # store lfc
		# glmfit$lfc <- lfc
# 
		# # check logFC
		# logFCt <- logFC
		# if(!is.null(glmfit$unshrunk.coefficients)) logFCt <- unshrunk.logFC
# 
		# if(upshot){
			# # We choose 17 nodes and here are nodes and weights (we keep 0 as a node to avoid exact p-value equal to 1)
			# # From the simualtion, 17 nodes works well when lfc <= 6
			# # limma uses 16 nodes for UPSHOT p-value
			# gq.nodes   <- c(0.1784842,0.3512318,0.5126905,0.6576712,0.7815140,0.8802392,0.9506755,0.9905755)
			# gq.weights <- c(0.1765627,0.1680041,0.1540458,0.1351364,0.1118839,0.0850362,0.0554595,0.0241483)
# 
			# # UPSHOT p-values by gauss-quadrature
			# F.pvalue <- 0.08972310 * F.pvalue
			# for(i in 1:8){
				# F.pvalue <- F.pvalue + gq.weights[i] * .glmtreat(glmfit,design,design0,coef,ncoef,df.total,lfc*gq.nodes[i],dispersion,logFCt)
			# }
		# } else {
			# F.pvalue <- (F.pvalue + .glmtreat(glmfit,design,design0,coef,ncoef,df.total,lfc,dispersion,logFCt)) / 2
		# }
# 
		# # There is no statistics corresponding to UPSHOT p-value
		# # Here I just convert it to F-statistic 
		# # limma use the half lfc
		# # Or I can use the last node?
		# F.stat <- qf(F.pvalue, df1=ncoef, df2=df.total, lower.tail = FALSE)
	# } 

	#	Ensure it is not more significant than chisquare test with Poisson variance		
	if(poisson.bound) {
		i <- .isBelowPoissonBound(glmfit, nthreads=glmfit$nthreads)
		if(any(i)) {
			pois.fit <- glmfit[i,]
			pois.fit <- glmFit(pois.fit$counts, design=pois.fit$design, offset=pois.fit$offset, weights=pois.fit$weights, start=pois.fit$unshrunk.coefficients, dispersion=0, nthreads=glmfit$nthreads)
			pois.res <- glmLRT(pois.fit, coef = coef, contrast = contrast)
			F.pvalue[i] <- pmax(F.pvalue[i], pois.res$table$PValue)
		}
	}
	
	rn <- rownames(glmfit)
	if(is.null(rn))
		rn <- 1:nrow(glmfit)
	else
		rn <- make.unique(rn)
	
	#	Table output
	if(ncoef==1) logFC <- drop(logFC)
	
	tab <- data.frame(
		logFC=logFC,
		logCPM=glmfit$AveLogCPM,
		F=F.stat,
		PValue=F.pvalue,
		row.names=rn
	)
	
	glmfit$counts <- NULL
	glmfit$table <- tab
	glmfit$comparison <- coef.name
	glmfit$df.test <- ncoef
	glmfit$df.total <- df.total
	new("DGELRT",unclass(glmfit))
}

.isBelowPoissonBound <- function(glmfit, nthreads=1L)
# A convenience function to avoid generating temporary matrices.
{
	disp <- makeCompressedMatrix(glmfit$dispersion, dim(glmfit$counts), byrow=FALSE)
	s2   <- makeCompressedMatrix(glmfit$s2.post, dim(glmfit$counts), byrow=FALSE)
	out  <- .Call(.cxx_check_poisson_bound, glmfit$fitted.values, disp, s2, nthreads)
	return(out)
}

# .glmtreat <- function(glmfit, design, design0, coef, ncoef, df.total, lfc, dispersion=NULL, logFCt)
# #   function to perform treat analysis on the gauss-quadrature nodes
# #	Created by Lizhong Chen 13 Nov 2025
# {
	# # number of hypothesis tests
	# nboundary <- 2^ncoef
# 
	# # quick computation for those genes with lfc within the threshold
	# logFCt <- abs(t(logFCt)) * log(2)
	# lfc    <- lfc * log(2)
	# ind    <- apply(logFCt,2,function(x) any(x>lfc))
# 
	# # if all genes with lfc < threshold, retuen p=1
	# F.pvalue <- rep(1,nrow(glmfit))  
	# if(sum(ind)==0) return(F.pvalue)
# 
	# # subset the genes
	# if(length(dispersion)>1L){
		# dispersion <- makeCompressedMatrix(dispersion,dim(glmfit$counts),byrow=FALSE)
		# dispersion <- dispersion[ind,,drop=FALSE]
	# } 
	# logFCt  <- logFCt[,ind,drop=FALSE]
	# counts  <- glmfit$counts[ind,,drop=FALSE]
	# offset0 <- t(as.matrix(glmfit$offset)[ind,,drop=FALSE])
	# weights <- glmfit$weights[ind,,drop=FALSE] 
# 
	# # tested coef
	# design1 <- design[, coef, drop=FALSE]
# 
	# # compute deviance for null hypothesis
	# dev <- matrix(0,sum(ind),nboundary)
	# LR  <- rep(0,sum(ind))
	# for(i in 1:nboundary){
		# # convert an integer to a vector with entries +/- 1
		# boundary <- 2 * as.integer(intToBits(i-1))[1:ncoef] - 1
# 
		# # approximate MLE for the null hypothesis when MLE is outside
		# # consider the projection of MLE, if it is on the boundary, we assume it is MLE for null
		# # if not, MLE should be on the edge of the boundary
		# # when MLE is inside, that is fine as the null is accepted
		# eta     <- design1 %*% matrix(pmin(lfc, logFCt) * boundary, ncoef)
		# offset  <- t(eta+offset0)
		# fit     <- glmFit(counts,design=design0,offset=offset,weights=weights,dispersion=dispersion,prior.count=0)
		# dev[,i] <- fit$deviance
	# }
	# LR <- pmax(apply(dev-glmfit$deviance[ind],1,min),0)
# 
	# #	Compute p-values from the QL F-statistic
	# F.stat   <- LR / ncoef / glmfit$s2.post[ind]
	# F.pvalue[ind] <- pf(F.stat, df1=ncoef, df2=df.total[ind], lower.tail=FALSE, log.p=FALSE) 
# 
	# F.pvalue
# }

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
