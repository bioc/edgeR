diffSpliceDGE <- function(glmfit, coef=ncol(glmfit$design), contrast=NULL, geneid, exonid=NULL, prior.count=0.125, robust=NULL, verbose=TRUE)
# Identify exons and genes with splice variants using negative binomial GLMs
# Yunshun Chen, Lizhong Chen and Gordon Smyth
# Created 29 March 2014.  Last modified 5 May 2024.
{
#	Check if glmfit is from glmFit() or glmQLFit()
	isLRT <- is.null(glmfit$df.prior)

#	keep robust consistent with glmfit object
	if(is.null(robust) & !isLRT){
		robust <- length(glmfit$df.prior) > 1L
	}

#	Check input (from diffSplice in limma)
	exon.genes <- glmfit$genes
	nexons <- nrow(glmfit)
	design <- glmfit$design

	if(is.null(exon.genes)) exon.genes <- data.frame(ExonID=1:nrow(glmfit))
	if(length(geneid)==1) {
		genecolname <- as.character(geneid)
		geneid <- exon.genes[[genecolname]]
	} else {
		exon.genes$GeneID <- geneid
		genecolname <- "GeneID"
	}
	if(!is.null(exonid))
		if(length(exonid)==1) {
			exoncolname <- as.character(exonid)
			exonid <- exon.genes[[exoncolname]]
		} else {
			exon.genes$ExonID <- exonid
			exoncolname <- "ExonID"
		}
	else
		exoncolname <- NULL

#	Sort by geneid
	if(is.null(exonid))
		o <- order(geneid)
	else
		o <- order(geneid,exonid)
	geneid <- geneid[o]
	exon.genes <- exon.genes[o,,drop=FALSE]
	glmfit <- glmfit[o, ]

#	Check design matrix
	design <- as.matrix(glmfit$design)
	nbeta <- ncol(design)
	if(nbeta < 2) stop("Need at least two columns for design, usually the first is the intercept column")
	coef.names <- colnames(design)

	coefficients <- glmfit$coefficients

#	Evaluate beta to be tested
#	Note that contrast takes precedence over coef: if contrast is given
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
		else
			coef.name <- coef.names[coef]
		beta <- coefficients[, coef, drop=FALSE]
#		g.ind <- rowSums(design[,-coef,drop=FALSE]) == 1L
	} else {
		contrast <- as.matrix(contrast)
		if(ncol(contrast) > 1L) {
			warning("contrast is a matrix, should be a vector. Using first column only.")
			contrast <- contrast[,1,drop=FALSE]
		}
		reform <- contrastAsCoef(design, contrast=contrast, first=TRUE)
		coef <- 1
		beta <- drop(coefficients %*% contrast)
		contrast <- drop(contrast)
		i <- contrast!=0
		coef.name <- paste(paste(contrast[i],coef.names[i],sep="*"),collapse=" ")
		design <- reform$design
#		g.ind <- abs(design[,1]) > 0
	}
	beta <- as.vector(beta)	
	
#	Null design matrix
	design0 <- design[, -coef, drop=FALSE]

# 	Count exons and get genewise variances
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
	ngenes <- sum(gene.keep)
	if(ngenes==0) stop("No genes with more than one exon")

	exon.keep   <- rep(gene.keep, gene.nexons)
	geneid      <- geneid[exon.keep]
	exon.genes  <- exon.genes[exon.keep, , drop=FALSE]
	beta        <- beta[exon.keep]
	gene.nexons <- gene.nexons[gene.keep]
	
#	Gene level information
	g <- rep(1:ngenes, times=gene.nexons)
	glmfit <- glmfit[exon.keep, ]
	gene.counts  <- rowsum(glmfit$counts, geneid, reorder=FALSE)
	fit.gene     <- glmFit(gene.counts, design, dispersion=0.05, offset=as.vector(glmfit$offset[1,]), prior.count=prior.count)
	gene.betabar <- fit.gene$coefficients[g, coef, drop=FALSE]

#	New offset
	offset.new   <- makeCompressedMatrix(glmfit$offset, dim(glmfit$counts), byrow=TRUE) + 
                    makeCompressedMatrix(tcrossprod(gene.betabar, design[,coef,drop=FALSE]))
	coefficients <- beta - gene.betabar

#	adjust dispersion for new QL method scaled by average QL dispersion
	if(is.null(glmfit$average.ql.dispersion)){
		dispersion <- glmfit$dispersion
	} else {
		dispersion <- glmfit$dispersion/glmfit$average.ql.dispersion
	}

#	Testing: Null model
	fit0 <- glmFit(glmfit$counts, design=design0, offset=offset.new, dispersion=dispersion, weights=glmfit$weights)

#	Prepare statistics summary
	exon.LR <- fit0$deviance - glmfit$deviance
	gene.LR <- rowsum(exon.LR, geneid, reorder=FALSE)
	exon.df.test <- fit0$df.residual - glmfit$df.residual
	gene.df.test <- rowsum(exon.df.test, geneid, reorder=FALSE) - 1

#	degree of freedom and deviance for QL methods
	if(!is.null(glmfit$df.residual.zeros)){
		exon.df.residual <- glmfit$df.residual.zeros
		exon.deviance <- glmfit$deviance
	}

	if(!is.null(glmfit$df.residual.adj)){
		exon.df.residual <- glmfit$df.residual.adj
		exon.deviance <- glmfit$deviance.adj
	}

#	exon and gene level test
	if(isLRT){
		exon.p.value <- pchisq(exon.LR, df=exon.df.test, lower.tail=FALSE, log.p=FALSE)
		gene.p.value <- pchisq(gene.LR, df=gene.df.test, lower.tail=FALSE, log.p=FALSE)
	} else {
		gene.df.residual <- rowsum(exon.df.residual, geneid, reorder=FALSE)
		gene.s2 <- rowsum(exon.deviance, geneid, reorder=FALSE) / gene.df.residual
		gene.df.residual[gene.df.residual < 0.99] <- 0

		squeeze <- squeezeVar(var=gene.s2, df=gene.df.residual, robust=robust)
		gene.df.total <- gene.df.residual + squeeze$df.prior
		gene.df.total <- pmin(gene.df.total, sum(gene.df.residual))
		gene.s2.post  <- squeeze$var.post

		exon.F <- exon.LR / exon.df.test / gene.s2.post[g]
		gene.F <- rowsum(exon.LR, geneid, reorder=FALSE) / gene.df.test / gene.s2.post
		exon.p.value <- pf(exon.F, df1=exon.df.test, df2=gene.df.total[g], lower.tail=FALSE, log.p=FALSE)
		gene.p.value <- pf(gene.F, df1=gene.df.test, df2=gene.df.total, lower.tail=FALSE, log.p=FALSE)

#		Ensure is not more significant than chisquare test for previous QL method
		i <- gene.s2.post[g] < 1
		if(any(i) & !is.null(glmfit$df.residual.zeros)) {
			chisq.pvalue    <- pchisq(exon.LR[i], df=exon.df.test[i], lower.tail=FALSE, log.p=FALSE)
			exon.p.value[i] <- pmax(exon.p.value[i], chisq.pvalue)
		}
	}

#	Gene Simes' p-values
#	exon.zeros <- rowSums(glmfit$counts[,g.ind]) > 0
	o <- order(g, exon.p.value, decreasing=FALSE)
	p <- exon.p.value[o]
	q <- rep(1, sum(gene.nexons))
	r <- cumsum(q) - rep(cumsum(q)[cumsum(gene.nexons)]-gene.nexons, gene.nexons)
#	pp <- p*rep(gene.nexons, gene.nexons)/r
	pp <- p*pmax(rep(gene.nexons-1, gene.nexons)/r,1)
	oo <- order(-g, pmin(pp,1), decreasing=TRUE)
	gene.Simes.p.value <- pp[oo][cumsum(gene.nexons)]

#	Output
	out <- new("DGELRT",list())
	out$comparison <- colnames(design)[coef]
	out$design <- design
	out$coefficients <- as.vector(coefficients)
	out$genes <- exon.genes
	out$genecolname <- genecolname
	out$exoncolname <- exoncolname
	
#	Exon level output
	out$exon.df.test <- exon.df.test
	if(isLRT){
		out$exon.LR <- exon.LR
	} else {
		out$exon.F <- exon.F
	}
	out$exon.p.value <- exon.p.value

#	Gene level output
	out$gene.df.test <- gene.df.test
	if(isLRT){
		out$gene.LR <- gene.LR
	} else {
		out$gene.df.prior <- squeeze$df.prior
		out$gene.df.residual <- gene.df.residual
		out$gene.F <- gene.F
	}
	out$gene.p.value <- gene.p.value
	out$gene.Simes.p.value <- gene.Simes.p.value

#	Which columns of exon.genes contain gene level annotation? (from diffSplice in limma)
	exon.lastexon <- cumsum(gene.nexons)
	exon.firstexon <- exon.lastexon-gene.nexons+1
	no <- logical(nrow(exon.genes))
	isdup <- vapply(exon.genes,duplicated,no)[-exon.firstexon,,drop=FALSE]
	isgenelevel <- apply(isdup,2,all)
	out$gene.genes <- exon.genes[exon.lastexon,isgenelevel, drop=FALSE]
	out$gene.genes$NExons <- gene.nexons

	out
}


topSpliceDGE <- function(lrt, test="Simes", number=10, FDR=1)
# Yunshun Chen and Gordon Smyth
# Created 29 March 2014.  Last modified 25 September 2015. 
{
	test <- match.arg(test,c("Simes","simes","gene","exon"))
	if(test=="simes") test <- "Simes"
	if(test=="exon") {
		number <- min(number, nrow(lrt$genes))
		P <- lrt$exon.p.value
		BH <- p.adjust(P, method="BH")
		if(FDR<1) number <- min(number, sum(BH<FDR))
		o <- order(P)[1:number]
		if(is.null(lrt$exon.F)){
			data.frame(lrt$genes[o,,drop=FALSE],logFC=lrt$coefficients[o],exon.LR=lrt$exon.LR[o],P.Value=P[o],FDR=BH[o])
		} else {
			data.frame(lrt$genes[o,,drop=FALSE],logFC=lrt$coefficients[o],exon.F=lrt$exon.F[o],P.Value=P[o],FDR=BH[o])
		}
	} else {
		number <- min(number, nrow(lrt$gene.genes))
		if(test=="Simes") P <- lrt$gene.Simes.p.value else P <- lrt$gene.p.value 
		BH <- p.adjust(P, method="BH")
		if(FDR<1) number <- min(number,sum(BH<FDR))
		o <- order(P)[1:number]
		if(test=="Simes"){
			data.frame(lrt$gene.genes[o,,drop=FALSE],P.Value=P[o],FDR=BH[o])
		} else {
			if(is.null(lrt$gene.F)){
				data.frame(lrt$gene.genes[o,,drop=FALSE],gene.LR=lrt$gene.LR[o],P.Value=P[o],FDR=BH[o])
			} else {
				data.frame(lrt$gene.genes[o,,drop=FALSE],gene.F=lrt$gene.F[o],P.Value=P[o],FDR=BH[o])
			}
		}
	}
}


plotSpliceDGE <- function(lrt, geneid=NULL, genecolname=NULL, rank=1L, FDR = 0.05)
# Plot exons of most differentially spliced gene
# Yunshun Chen and Gordon Smyth
# Created 29 March 2014.  Last modified 17 April 2024.
{
	if(is.null(genecolname)) 
		genecolname <- lrt$genecolname
	else
		genecolname <- as.character(genecolname)
	
	if(is.null(geneid)) {
		if(rank==1L)
			i <- which.min(lrt$gene.Simes.p.value)
		else
			i <- order(lrt$gene.Simes.p.value)[rank]
		geneid <- paste(lrt$gene.genes[i, genecolname], collapse = ".")
	} else {
		geneid <- as.character(geneid)
		i <- which(lrt$gene.genes[, genecolname]==geneid)[1]
		if(!length(i)) stop(paste("geneid",geneid,"not found"))
	}

	exon.lastexon <- cumsum(lrt$gene.genes$NExons[1:i])
	j <- (exon.lastexon[i]-lrt$gene.genes$NExons[i]+1):exon.lastexon[i]

	exoncolname <- lrt$exoncolname
	if(is.null(exoncolname)){
		plot(lrt$coefficients[j], xlab="Exon", ylab="logFC (this exon vs the average)", main=geneid, type="b")
	}
	# Plot exons and mark exon ids on the axis
	if(!is.null(exoncolname)) {
		exon.id <- lrt$genes[j, exoncolname]
		xlab <- paste("Exon", exoncolname, sep=" ")
		
		plot(lrt$coefficients[j], xlab="", ylab="logFC (this exon vs the average)", main=geneid, type="b", xaxt="n")
		axis(1, at=1:length(j), labels=exon.id, las=2, cex.axis=0.6)
		mtext(xlab, side=1, padj=5.2)

		# Mark the topSpliced exons
		top <- topSpliceDGE(lrt, number=Inf, test="exon", FDR=FDR)
		m <- which(top[,genecolname] %in% lrt$gene.genes[i,genecolname])

		if(length(m) > 0){
			if(length(m) == 1) cex <- 1.5 else{
				abs.fdr <- abs(log10(top$FDR[m]))
				from <- range(abs.fdr)
				to <- c(1,2)
				if(diff(from) == 0)
					cex <- rep(1.5, length(m))
				else
					cex <- (abs.fdr - from[1])/diff(from) * diff(to) + to[1]
			}	
			mark <- match(top[m, exoncolname], exon.id)
			points((1:length(j))[mark], lrt$coefficients[j[mark]], col = "red", pch = 16, cex = cex)
		}
	}
	abline(h=0,lty=2)
	invisible()
}
