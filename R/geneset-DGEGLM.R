roast.DGEGLM <- function(y, index=NULL, design=NULL, contrast=ncol(design), geneid=NULL, set.statistic="mean", gene.weights=NULL, nrot=1999, ...)
#	Rotation gene set testing for a DGEGLM class object
#	Yunshun Chen
#	Created 12 Oct 2023.
{
	if(is.null(design)) design <- y$design
	if(length(geneid)==1L) geneid <- y$genes[[geneid]]
	y <- .zscoreGLM(y=y, design=design, contrast=contrast)
	roast(y=y, index=index, design=design, contrast=contrast, geneid=geneid, set.statistic=set.statistic, gene.weights=gene.weights, var.prior=1, df.prior=Inf, nrot=nrot, approx.zscore=TRUE)
}

mroast.DGEGLM <- function(y, index=NULL, design=NULL, contrast=ncol(design), geneid=NULL, set.statistic="mean", gene.weights=NULL, nrot=1999, adjust.method="BH", midp=TRUE, sort="directional", ...)
#	Rotation gene set testing for a DGEGLM class object (multiple sets)
#	Yunshun Chen
#	Created 12 Oct 2023.
{
	if(is.null(design)) design <- y$design
	if(length(geneid)==1L) geneid <- y$genes[[geneid]]
	y <- .zscoreGLM(y=y, design=design, contrast=contrast)
	mroast(y=y, index=index, design=design, contrast=contrast, geneid=geneid, set.statistic=set.statistic, gene.weights=gene.weights, var.prior=1, df.prior=Inf, nrot=nrot, approx.zscore=TRUE, adjust.method=adjust.method, midp=midp, sort=sort)
}

fry.DGEGLM <- function(y, index=NULL, design=NULL, contrast=ncol(design), geneid=NULL, sort="directional", ...)
#	Rotation gene set testing for a DGEGLM class object (multiple sets)
#	Yunshun Chen
#	Created 12 Oct 2023
{
	if(is.null(design)) design <- y$design
	if(length(geneid)==1L) geneid <- y$genes[[geneid]]
	y <- .zscoreGLM(y=y, design=design, contrast=contrast)
	fry(y=y, index=index, design=design, contrast=contrast, standardize="none", geneid=geneid, sort=sort)
}

camera.DGEGLM <- function(y, index, design=NULL, contrast=ncol(design), weights=NULL, use.ranks=FALSE, allow.neg.cor=FALSE, inter.gene.cor=0.01, sort=TRUE, ...)
#	Rotation gene set testing for a DGEGLM class object accounting for inter-gene correlation
#	Yunshun Chen
#	Created 12 Oct 2023.
{
	if(is.null(design)) design <- y$design
	y <- .zscoreGLM(y=y, design=design, contrast=contrast)
	camera(y=y, index=index, design=design, contrast=contrast, weights=weights, use.ranks=use.ranks, allow.neg.cor=allow.neg.cor, inter.gene.cor=inter.gene.cor, trend.var=FALSE, sort=sort)
}

romer.DGEGLM <- function(y, index, design=NULL, contrast=ncol(design), ...)
#	rotation mean-rank version of GSEA (gene set enrichment analysis) for a DGEGLM class object
#	Yunshun Chen
#	Created 12 Oct 2023.
{
	if(is.null(design)) design <- y$design
	y <- .zscoreGLM(y=y, design=design, contrast=contrast)
	romer(y=y, index=index, design=design, contrast=contrast, ...)
}


.zscoreGLM <- function(y, design=NULL, contrast=ncol(design))
#	Calculate NB z-scores given a DGEGLM object.
#	Yunshun Chen and Lizhong Chen
#	Created 19 Oct 2023.
{
#	Check for all zero counts
	allzero <- rowSums(y$counts>1e-8)==0
	if(any(allzero)) warning(sum(allzero),"rows with all zero counts")

#	Prepare count matrix
	counts <- y$counts

#	check average QL dispersion
	if(!is.null(y$average.ql.dispersion)){
		counts <- counts / pmax(1, y$var.prior)
	}

#	Make default design matrix from group factor
	if(is.null(design)) design <- y$design
	nbeta <- ncol(design)
	if(nbeta < 2) stop("design matrix must have at least two columns")

#	contrast could be a coef name
	if(is.character(contrast)) {
		if(length(contrast)>1) stop("contrast should specify only one column of design")
		contrast <- which(contrast==colnames(design))
		if(!length(contrast)) stop("contrast doesn't match any column of design")
	}

#	Construct null hypothesis design matrix
	if(length(contrast) == 1) {
		design0 <- design[,-contrast,drop=FALSE]
	} else {
		design <- contrastAsCoef(design,contrast=contrast,first=FALSE)$design
		design0 <- design[,-nbeta,drop=FALSE]
	}

#	Null hypothesis fit
#	the scale counts for new QL method still follow the original dispersion
	fit.null <- glmFit(counts, design=design0, dispersion=y$dispersion, offset=y$offset, weights=y$weights, prior.count=0)

#	Quantile residuals from null fit
#	Applying a floor to mu avoids infinite z-scores when mu=0
	y <- zscoreNBinom(counts, mu=pmax(fit.null$fitted.values,1e-17), size=1/y$dispersion)
	y
}
