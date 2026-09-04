#	Fit Binomial Models with adjusted deviance statistics
binQLFit <- function(y, ...)
	UseMethod("binQLFit")

binQLFit.PCList <- function(y, design=NULL, abundance.trend=TRUE, covariate.trend=NULL, robust=TRUE, ...)
#	Created 11 Jul 2026. Last modified 11 Jul 2026.
{
	if(is.null(design)) {
		design <- y$design
		if(is.null(design)) {
			group <- droplevels(as.factor(y$samples$group))
			if(nlevels(group) > 1L) design <- model.matrix(~y$samples$group)
		}
	}

	if(is.null(y$AveLogCPM)) y$AveLogCPM <- aveLogCPM(y)
	
	fit <- binQLFit.default(y=y$counts,z=y$counts2,design=design,weights=y$weights,AveLogCPM=y$AveLogCPM,abundance.trend=abundance.trend,covariate.trend=covariate.trend,robust=robust, ...)
	
	fit$samples   <- y$samples
	fit$genes     <- y$genes
	fit$AveLogCPM <- y$AveLogCPM
	fit
}

binQLFit.default <- function(y, z, design=NULL, weights=NULL, offset=NULL, AveLogCPM=NULL, abundance.trend=TRUE, covariate.trend=NULL, robust=TRUE, prior.count=1, nthreads=1L, ...)
#	Fit binomial generalized linear model for each feature
#	Lizhong Chen, Gordon Smyth
#	Created 11 Jul 2026. Last modified 11 Jul 2026.
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

	#	compressed offset (do not overwrite 'offset'; it is reused by the prior.count recursion below)
	if(is.null(offset)) offset0 <- makeCompressedMatrix(0, dim(y)) else offset0 <- .compressOffsets(y, offset)

	#	Fit + QL adjustment (leverage-adjusted deviance/df/quasi-dispersion) in a
	#	SINGLE C call, replacing .cxx_bin_fit + .cxx_compute_adj_vec_bin.  Coverage,
	#	proportion and coverage-scaled weights are formed inside C from the raw counts.
	fit <- .Call(.cxx_bin_ql_fit, y, z, weights, offset0, design, 50L, 1e-8, 100L, 1e-6, nthreads)
	s2  <- fit$s2
	df.adj <- fit$df.residual.adj
	fit$s2 <- NULL

	# --- replaced two-call path (kept for reference) ---
	# coverage <- y + z
	# weights0 <- if(is.null(weights)) coverage else as.matrix(weights) * coverage
	# prop <- y / pmax(coverage,1)
	# fit <- .Call(.cxx_bin_fit, prop, offset0, weights0, design, 50L, 1e-8, 100L, 1e-6, nthreads)
	# out <- .Call(.cxx_compute_adj_vec_bin, prop, fit$fitted.values, design, coverage, weights0, nthreads)
	# s2  <- out$s2
	# fit$df.residual.adj <- df.adj <- out$df
	# fit$deviance.adj <- out$deviance

	dimnames(fit$coefficients)  <- list(rownames(y), colnames(design))
	dimnames(fit$fitted.values) <- dimnames(y)
#	Average log2-CPM of the total coverage (as in aveLogCPM.PCList)
	if(is.null(AveLogCPM))
		AveLogCPM <- aveLogCPM(y+z, lib.size=colSums(y)+colSums(z), prior.count=2, dispersion=0, weights=weights)
	fit$AveLogCPM <- AveLogCPM

#	Covariate for trended prior for quasi-dispersion
	if(is.null(covariate.trend)) {
		if(!abundance.trend) AveLogCPM <- NULL
	} else {
		AveLogCPM <- covariate.trend
	}
	
#	Empirical Bayes moderation of quasi-likelihood dispersions
	s2.fit <- squeezeVar(s2,df=df.adj,covariate=AveLogCPM,robust=robust,legacy=FALSE)
	
#	Storing results
	fit$df.prior  <- s2.fit$df.prior
	fit$s2.post   <- s2.fit$var.post
	fit$s2.prior  <- s2.fit$var.prior  
	fit$average.ql.dispersion <- median(s2.fit$var.prior)

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

binQLFTest <- function(glmfit, coef=ncol(glmfit$design), contrast=NULL, or = 1, lor = NULL, upshot=TRUE)
#	Quasi-likelihood F-tests for binomial models.
#	Lizhong Chen, Gordon Smyth
#	Created 09 Jul 2026. Last modified 11 Jul 2026
{
	if(!is(glmfit,"DGEBIN")) stop("glmfit must be a DGEBIN object produced by binQLFit")
	if(is.null(glmfit$s2.post)) stop("need to run binQLFit before binQLFTest")
	nlibs <- ncol(glmfit)

	#	Check design matrix
	design <- as.matrix(glmfit$design)
	nbeta  <- ncol(design)
	if(nbeta < 2) stop("Need at least two columns for design, usually the first is the intercept column")
	coef.names <- colnames(design)
	
	#	Evaluate logOR for coef to be tested
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
		logOR <- glmfit$coefficients[,coef,drop=FALSE]/log(2)
		if(!is.null(glmfit$unshrunk.coefficients)) unshrunk.logOR <- glmfit$unshrunk.coefficients[,coef,drop=FALSE]/log(2)
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
		logOR  <- (glmfit$coefficients %*% contrast)/log(2)
		if(!is.null(glmfit$unshrunk.coefficients)) unshrunk.logOR <- (glmfit$unshrunk.coefficients %*% contrast)/log(2)
	}
	
	# number of coefficients
	ncoef <- length(coef)

	# Null design matrix
	design0 <- design[, -coef, drop=FALSE]

	# check lor 
	if(is.null(lor))    lor <- log2(or)
	if(length(lor)==1L) lor <- rep(lor, ncoef)
	if(sum(lor < 0) > 0) stop("lor should be non-negative")
	if((length(lor)!=ncoef) & (length(lor)>1)) stop("'or' or 'lor' vector of wrong length")

	# starting p-values for normal F-test    
	# fit <- binFit(glmfit$counts,glmfit$counts2,design=design0,weights=glmfit$weights,bias.correction=FALSE,prior.count=0)
	# LR  <- pmax(fit$deviance-glmfit$deviance,0)

	# degree of freedom
	df.residual <- glmfit$df.residual.adj
	df.total    <- glmfit$df.prior + df.residual
	df.residual.total <- sum(glmfit$df.residual)
	df.total <- pmin(df.total, df.residual.total)

	#	QL F-statistic + p-value + UPSHOT/TREAT computation in a SINGLE C call.
	#	The base null refit, F-statistic/pf, the TREAT Gauss-quadrature (formerly
	#	the R helper .bintreat) and the final qf all run in C via .cxx_bin_ftest,
	#	reusing the same binomial fitter (bin_fit_mat) that binFit drives through
	#	.cxx_bin_fit, so results match the former R path to ~1e-8.
	#	The logORt selection and the glmfit$lor side-effect stay in R; nthreads is
	#	1L because the former R TREAT refits (via binFit) ran single-threaded.
	logORt <- logOR
	if(!is.null(glmfit$unshrunk.coefficients)) logORt <- unshrunk.logOR
	out <- .Call(.cxx_bin_ftest, glmfit$counts, glmfit$counts2, glmfit$weights, design, coef, glmfit$deviance, glmfit$s2.post, df.total, lor, logORt, upshot, glmfit$nthreads)
	F.stat   <- out$F
	F.pvalue <- out$PValue

	#--- ORIGINAL F-statistic + TREAT/UPSHOT + qf (replaced by .cxx_bin_ftest above; kept for reference) ---
	# #	Compute p-values from the QL F-statistic
	# F.stat   <- LR / ncoef / glmfit$s2.post
	# F.pvalue <- pf(F.stat, df1=ncoef, df2=df.total, lower.tail=FALSE, log.p=FALSE)
	#
	# # perform treat analysis if any lor > 0
	# if(any(lor>0)){
	# 	# store lor
	# 	glmfit$lor <- lor
	#
	# 	# check logOR
	# 	logORt <- logOR
	# 	if(!is.null(glmfit$unshrunk.coefficients)) logORt <- unshrunk.logOR
	#
	# 	if(upshot){
	# 		# We choose 17 nodes and here are nodes and weights (we keep 0 as a node to avoid exact p-value equal to 1)
	# 		# From the simualtion, 17 nodes works well when lor <= 6
	# 		# limma uses 16 nodes for UPSHOT p-value
	# 		gq.nodes   <- c(0.1784842,0.3512318,0.5126905,0.6576712,0.7815140,0.8802392,0.9506755,0.9905755)
	# 		gq.weights <- c(0.1765627,0.1680041,0.1540458,0.1351364,0.1118839,0.0850362,0.0554595,0.0241483)
	#
	# 		# UPSHOT p-values by gauss-quadrature
	# 		F.pvalue <- 0.08972310 * F.pvalue
	# 		for(i in 1:8){
	# 			F.pvalue <- F.pvalue + gq.weights[i] * .bintreat(glmfit,design,design0,coef,ncoef,df.total,lor*gq.nodes[i],logORt)
	# 		}
	# 	} else {
	# 		F.pvalue <- (F.pvalue + .bintreat(glmfit,design,design0,coef,ncoef,df.total,lor,logORt)) / 2
	# 	}
	#
	# 	# There is no statistics corresponding to UPSHOT p-value
	# 	# Here I just convert it to F-statistic
	# 	# limma use the half lor
	# 	# Or I can use the last node?
	# 	F.stat <- qf(F.pvalue, df1=ncoef, df2=df.total, lower.tail = FALSE)
	# }
	
	rn <- rownames(glmfit)
	if(is.null(rn))
		rn <- 1:nrow(glmfit)
	else
		rn <- make.unique(rn)
	
	#	Table output
	if(ncoef==1) logOR <- drop(logOR)
	
	tab <- data.frame(
		logOR=logOR,
		logCPM=glmfit$AveLogCPM,
		F=F.stat,
		PValue=F.pvalue,
		row.names=rn
	)
	
	if(any(lor>0)) glmfit$lor <- lor
	glmfit$counts <- NULL
	glmfit$counts2 <- NULL
	glmfit$table <- tab
	glmfit$comparison <- coef.name
	glmfit$df.test <- ncoef
	glmfit$df.total <- df.total
	new("DGELRT",unclass(glmfit))
}

# .bintreat <- function(glmfit, design, design0, coef, ncoef, df.total, lor, logORt)
# #   function to perform treat analysis on the gauss-quadrature nodes
# #	Created by Lizhong Chen 13 Nov 2025
# {
	# # number of hypothesis tests
	# nboundary <- 2^ncoef
# 
	# # quick computation for those genes with lor within the threshold
	# logORt <- abs(t(logORt)) * log(2)
	# lor    <- lor * log(2)
	# ind    <- apply(logORt,2,function(x) any(x>lor))
# 
	# # if all genes with lor < threshold, retuen p=1
	# F.pvalue <- rep(1,nrow(glmfit))  
	# if(sum(ind)==0) return(F.pvalue)
# 
	# # subset the genes
	# logORt  <- logORt[,ind,drop=FALSE]
	# counts  <- glmfit$counts[ind,,drop=FALSE]
	# counts2 <- glmfit$counts2[ind,,drop=FALSE]
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
		# eta    <- design1 %*% matrix(pmin(lor, logORt) * boundary, ncoef)
		# offset <- t(eta)
		# fit    <- binFit(counts,counts2,design=design0,offset=offset,weights=weights,bias.correction=FALSE,prior.count=0)
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