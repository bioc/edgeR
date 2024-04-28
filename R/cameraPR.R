cameraPR.DGELRT <- function(statistic, index, use.ranks=FALSE, inter.gene.cor=0.01, sort=TRUE, ...)
#	Camera gene set testing using pre-ranked results from a DGELRT class object
#	Gordono Smyth
#	Created 27 Apr 2024.
{
	if(statistic$df.test[1] > 1L) stop("camera testing is not implemented for tests on more than 1 df")
	if(!is.null(statistic$table$LR)) {
		z <- sign(statistic$table$logFC) * sqrt(statistic$table$LR)
	} else {
		if(is.null(statistic$table$F)) stop("Cannot find test statistic")
		tstat <- sign(statistic$table$logFC) * sqrt(statistic$table$F)
		z <- zscoreT(tstat,df=statistic$df.total,approx=TRUE)
	}
	cameraPR(z, index=index, use.ranks=use.ranks, inter.gene.cor=inter.gene.cor, sort=sort)
}
