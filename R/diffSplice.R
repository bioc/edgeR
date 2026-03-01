diffSplice.DGEGLM <- function(fit, coef=ncol(fit$design), contrast=NULL, geneid, exonid=NULL, robust=NULL, nexons.approx=10L, verbose=TRUE, ...)
# Identify exons and genes with splice variants using negative binomial GLMs
# Lizhong Chen, Yunshun Chen and Gordon Smyth
# Created 12 Feb 2025.  Last modified 29 Mar 2025.
{
	if(is.null(fit$s2.post)) stop("need to run glmQLFit before diffSplice")
	
#	By defautl, keep robust consistent with fit object
	if(is.null(robust)){
		robust <- length(fit$df.prior) > 1L
	}
	
#	Make sure there is always an annotation frame
	exon.genes <- fit$genes
	if(is.null(exon.genes)){
		exon.genes  <- data.frame(ExonID=1:nrow(fit))    
	} 
	
#	Get ID columns for genes and exons  
	if(length(geneid)==1) {
		genecolname <- as.character(geneid)
		geneid      <- exon.genes[[genecolname]]
	} else {
		exon.genes$GeneID <- geneid
		genecolname <- "GeneID"
	}
	if(!is.null(exonid)){
		if(length(exonid)==1) {
			exoncolname <- as.character(exonid)
			exonid <- exon.genes[[exoncolname]]
		} else {
			exon.genes$ExonID <- exonid
			exoncolname <- "ExonID"
		}
	}
	else{
		exoncolname <- NULL
	}
	
#	Treat NA geneids as genes with one exon
	if(anyNA(geneid)) {
		isna <- which(is.na(geneid))
		geneid[isna] <- paste0("NA",1:length(isna))
	}
	
#	Sort by geneid
	if(is.null(exonid)){
		o <- order(geneid)    
	} else {
		o <- order(geneid,exonid)
	}
	
	exon.genes <- exon.genes[o,,drop=FALSE]
	geneid <- geneid[o]
	fit <- fit[o, ]
	
#	Check design matrix
	design <- as.matrix(fit$design)
	nbeta  <- ncol(design)
	if(nbeta < 2) stop("Need at least two columns for design, usually the first is the intercept column")
	coef.names <- colnames(design)
	
#	Evaluate beta to be tested
#	Note that contrast takes precedence over coef, if contrast is given
#	then reform design matrix so that contrast of interest is the first column
	if(is.null(contrast)) {
		if(length(coef) > 1) {
			warning("coef is a vector, should be a single value. Using first value only.")
			coef <- coef[1]
		}
		if(is.character(coef)) {
			check.coef <- coef %in% colnames(design)
			if(any(!check.coef)) stop("One or more named coef arguments do not match a column of the design matrix.")
			coef.name <- coef
			coef <- match(coef, colnames(design))
		}
		else{
			coef.name <- coef.names[coef]
		}
		beta <- fit$coefficients
	} else {
		contrast <- as.matrix(contrast)
		if(ncol(contrast) > 1L) {
			warning("contrast is a matrix, should be a vector. Using first column only.")
			contrast <- contrast[,1,drop=FALSE]
		}
		reform <- contrastAsCoef(design, contrast=contrast, first=TRUE)
		coef   <- 1
		contrast  <- drop(contrast)
		i         <- contrast!=0
		coef.name <- paste(paste(contrast[i],coef.names[i],sep="*"),collapse=" ")
		design    <- reform$design
		refit     <- glmFit(fit$counts, design, dispersion=fit$dispersion, offset=fit$offset, weights=fit$weights)
		beta      <- refit$coefficients
	}
	
#	Count exons and get genewise variances
	nexons      <- nrow(fit)
	gene.nexons <- rowsum(rep(1,nexons), geneid, reorder=FALSE)
	if(verbose) {
		cat("Total number of exons: ", nexons, "\n")
		cat("Total number of genes: ", length(gene.nexons), "\n")
		cat("Number of genes with 1 exon: ", sum(gene.nexons==1), "\n")
		cat("Mean number of exons in a gene: ", round(mean(gene.nexons),0), "\n")
		cat("Max number of exons in a gene: ", max(gene.nexons), "\n")
	}
	
#	Remove genes with only 1 exon
	gene.keep <- gene.nexons > 1
	ngenes    <- sum(gene.keep)
	if(ngenes==0) stop("No genes with more than one exon")
	
	exon.keep   <- rep(gene.keep, gene.nexons)
	nexons      <- sum(exon.keep)
	
	geneid      <- geneid[exon.keep]
	exon.genes  <- exon.genes[exon.keep, , drop=FALSE]
	gene.nexons <- gene.nexons[gene.keep]
	exon.nexons <- rep(gene.nexons, times=gene.nexons)
	fit         <- fit[exon.keep, ]
	beta        <- beta[exon.keep, ]
	
#	Check weights
	weights <- .compressWeights(fit$counts,fit$weights)
	
#	Check offsets
	offset  <- .compressOffsets(fit$counts,fit$offset,lib.size=NULL)
	
#	Check dispersion, adjust dispersion for edgeR v4 QL method scaled by average QL dispersion
	if(is.null(fit$average.ql.dispersion)){
		dispersion <- fit$dispersion
	} else {
		dispersion <- fit$dispersion/fit$average.ql.dispersion
	}
	
#	Testing null model on gene level
	gene.lastexon  <- cumsum(gene.nexons)
	gene.firstexon <- gene.lastexon-gene.nexons+1
	names(gene.lastexon) <- names(gene.firstexon) <- geneid[gene.firstexon]

	gene.dev <- rowsum(fit$deviance, geneid, reorder=FALSE)
	exon.LR <- exon.coef <- matrix(0,nexons,1)
	gene.LR <- matrix(0,ngenes,1)
	colnames(exon.LR) <- colnames(exon.coef) <- colnames(gene.LR) <- coef.name
	rownames(exon.LR) <- rownames(exon.coef) <- exon.genes$ExonID
	rownames(gene.LR) <- geneid[gene.firstexon]

	for(nexon in unique(gene.nexons))
	{
	#	subset data according to number of transcripts
		exon.keep <- exon.nexons == nexon
		gene.keep <- gene.nexons == nexon
		ngene     <- sum(gene.keep)
		
		y0       <- fit$counts[exon.keep, ,drop=FALSE]
		offset0  <- offset[exon.keep, ,drop=FALSE]
		weights0 <- weights[exon.keep, ,drop=FALSE]
		beta0    <- beta[exon.keep, ,drop=FALSE]
		
		exon.dev <- fit$deviance[exon.keep]
		
		if(length(dispersion) > 1L){
			dispersion0 <- dispersion[exon.keep]
		} else {
			dispersion0 <- dispersion
		}
		
	#	fit on gene level
		fit0     <- .fitByGene(y0, design, coef, beta0, dispersion0, offset0, weights0, nexon)
		gene.LR[gene.keep] <- fit0$deviance - gene.dev[gene.keep]
		betabar  <- rep(fit0$beta, each=nexon)
		
	#	fit on transcript level
		if(nexon > nexons.approx){
			u0 <- matrix(t(fit0$fitted.values),ngene*nexon, ncol(y0), byrow = TRUE)
			exon.LR[exon.keep]   <- nbinomDeviance(y0,u0,dispersion=dispersion0, weights=weights0) - exon.dev
			exon.coef[exon.keep] <- beta0[,coef] - betabar
		} else if(nexon > 2L){
			i0 <- (0L:(ngene-1L))*nexon
			for(k in 1:nexon){
				jk <- i0+k
				y1 <- y0[-jk, ,drop=FALSE]
				offset1  <- offset0[-jk, ,drop=FALSE]
				weights1 <- weights0[-jk, ,drop=FALSE]
				beta1    <- beta0[-jk, ,drop=FALSE]
				if(length(dispersion0) > 1L){
					dispersion1 <- dispersion0[-jk]
				} else {
					dispersion1 <- dispersion0
				}       
				fit1 <- .fitByGene(y1, design, coef, beta1, dispersion1, offset1, weights1, nexon-1)
				exon.LR[exon.keep][jk]   <- fit0$deviance - (exon.dev[jk]+fit1$deviance)
				exon.coef[exon.keep][jk] <- beta0[jk,coef] - fit1$beta
			}
		} else {
			i0 <- (0L:(ngene-1L))*2 + 1
			exon.LR[exon.keep]          <- rep(gene.LR[gene.keep], each=2)
			exon.coef[exon.keep][i0]    <- beta0[i0, coef] - beta0[i0+1, coef]
			exon.coef[exon.keep][i0+1]  <- beta0[i0+1, coef] - beta0[i0, coef]
		}
}
	
#	Prepare statistics summary
	exon.df.test <- rep(1, nexons)
	gene.df.test <- gene.nexons - 1
	
	if(is.null(fit$average.ql.dispersion)){
		exon.stat <- cbind(fit$df.residual.zeros,fit$deviance)
	} else {
		exon.stat <- cbind(fit$df.residual.adj,fit$deviance.adj)
	}
	gene.sum  <- rowsum(exon.stat, geneid, reorder=FALSE)
	gene.df.residual <- gene.sum[,1]
	gene.s2          <- gene.sum[,2] / gene.sum[,1]
	
	squeeze       <- squeezeVar(var=gene.s2, df=gene.df.residual, robust=robust)
	gene.df.total <- gene.df.residual + squeeze$df.prior
	gene.df.total <- pmin(gene.df.total, sum(gene.df.residual))
	gene.s2.post  <- squeeze$var.post
	
	exon.df.total <- rep(gene.df.total, times=gene.nexons)
	exon.s2.post  <- rep(gene.s2.post,  times=gene.nexons)
	
	exon.F <- exon.LR / exon.df.test / exon.s2.post
	gene.F <- gene.LR / gene.df.test / gene.s2.post
	
	exon.t <- sqrt(pmax(0,exon.F)) * sign(exon.coef)

	exon.p.value <- pf(exon.F, df1=exon.df.test, df2=exon.df.total, lower.tail=FALSE, log.p=FALSE)
	gene.p.value <- pf(gene.F, df1=gene.df.test, df2=gene.df.total, lower.tail=FALSE, log.p=FALSE)
	
#	Output
	out <- new("MArrayLM",list())
	out$design       <- design
	out$comparison   <- colnames(design)[coef]
	out$genes        <- exon.genes
	out$genecolname  <- genecolname
	out$exoncolname  <- exoncolname
	out$coefficients <- exon.coef
	out$nexons.approx <- nexons.approx

#	Exon level output
	out$t       <- exon.t
	out$p.value <- exon.p.value
	
#	Gene level output
	out$gene.df.prior    <- squeeze$df.prior
	out$gene.df.residual <- gene.df.residual
	out$gene.df.total    <- gene.df.total
	out$gene.s2          <- gene.s2
	out$gene.s2.post     <- gene.s2.post
	out$gene.F             <- gene.F
	out$gene.F.p.value     <- gene.p.value
	
#	Which columns of exon.genes contain gene level annotation? (from diffSplice in limma)
	no <- logical(nrow(exon.genes))
	isdup <- vapply(exon.genes,duplicated,no)[-gene.firstexon,,drop=FALSE]
	isgenelevel <- apply(isdup,2,all)
	out$gene.genes <- exon.genes[gene.lastexon,isgenelevel, drop=FALSE]
	row.names(out$gene.genes) <- out$gene.genes[[genecolname]]
	out$gene.genes$NExons <- gene.nexons
	out$gene.firstexon <- gene.firstexon
	out$gene.lastexon  <- gene.lastexon
	
#	Gene Simes' p-values
	g <- rep(1:ngenes, times=gene.nexons)
	o <- order(g, exon.p.value, decreasing=FALSE)
	p <- exon.p.value[o]
	q <- rep(1, sum(gene.nexons))
	r <- cumsum(q) - rep(cumsum(q)[gene.lastexon]-gene.nexons, gene.nexons)
	pp <- p*rep(gene.nexons, gene.nexons)/r
	oo <- order(-g, pp, decreasing=TRUE)

	gene.simes.p.value <- gene.bonferroni.p.value <- gene.F
	
	gene.simes.p.value[,1] <- pp[oo][gene.lastexon]
	gene.bonferroni.p.value[,1] <- pmin(p[gene.firstexon]*(gene.nexons-1),1)
	
	out$gene.simes.p.value <- gene.simes.p.value
	out$gene.bonferroni.p.value <- gene.bonferroni.p.value
	
	out
}

.fitByGene <- function(counts, design, coef, beta, dispersion, offset, weights, transcripts.per.gene)
{
#	fit the null model for the genes with the same number of transcripts
#	Created by Gordon 11 Feb 2025
#	Modified by Lizhong Chen 12 Feb 2025

#	Number of genes and transcripts
	nsamples <- ncol(counts)
	ngenes   <- nrow(counts) / transcripts.per.gene
	
#	Design matrix for null model, first column to be tested
	design0 <- diag(transcripts.per.gene) %x% design[, -coef, drop=FALSE]
	design0 <- cbind(rep(design[,coef],transcripts.per.gene),design0)
	
#	Update counts, offsets, weights, starting values
	y       <- .convertMatGeneToTranscript(counts,  ngenes, nsamples*transcripts.per.gene, transcripts.per.gene)  
	offset  <- .convertMatGeneToTranscript(offset,  ngenes, nsamples*transcripts.per.gene, transcripts.per.gene)  
	weights <- .convertMatGeneToTranscript(weights, ngenes, nsamples*transcripts.per.gene, transcripts.per.gene) 
	beta    <- .convertBetaGeneToTranscript(beta, coef, ngenes, (ncol(beta)-1)*transcripts.per.gene, transcripts.per.gene)
	
#	Update dispersion
#	Note here we assume the dispersion is a scalar or a vector with length equal to number of rows
	if(length(dispersion) > 1L){
		dispersion <- matrix(rep(dispersion, each=nsamples), ngenes, nsamples*transcripts.per.gene, byrow = TRUE)
	}
	
#	fit null model
	fit <- mglmLevenberg(y, design0, dispersion=dispersion, offset=offset, weights=weights, coef.start=beta)
	
#	return deviance and average log fold change
	list(deviance = fit$deviance, beta = fit$coefficients[,1], fitted.values=fit$fitted.values)
}

.convertMatGeneToTranscript <- function(x, nrow, ncol, transcripts.per.gene)
{
#	convert a count matrix or CompressedMatrix on gene level into exon level
#	Created by Lizhong Chen, 12 Feb 2025
	out <- x
	if(inherits(x,"CompressedMatrix")){
		if(attr(x, "repeat.row") & attr(x, "repeat.col")){
			attr(out,"Dims") <- as.integer(c(nrow,ncol))
		} else if(attr(x, "repeat.row")){
			out <- makeCompressedMatrix(rep(as.vector(x), transcripts.per.gene), c(nrow, ncol))
		} else {
			out <- as.matrix(x)
			out <- matrix(t(out), nrow, ncol, byrow = TRUE)
		}
	} else {
		out <- as.matrix(x)
		out <- matrix(t(out), nrow, ncol, byrow = TRUE)    
	}
	out
}

.convertBetaGeneToTranscript <- function(beta, coef, nrow, ncol, transcripts.per.gene)
{
#	convert a coefficient matrix from gene to transcript as starting coefficients
#	Created by Lizhong Chen, 18 Feb 2025
	beta1 <- rowMeans(matrix(beta[,coef], nrow, transcripts.per.gene, byrow = TRUE))
	beta0 <- .convertMatGeneToTranscript(beta[,-coef, drop=FALSE], nrow, ncol, transcripts.per.gene)
	cbind(beta1, beta0)
}
