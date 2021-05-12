plotMDS.DGEList <- function (x,top=500,labels=NULL,pch=NULL,cex=1,dim.plot=c(1,2),gene.selection="pairwise",xlab=NULL,ylab=NULL,prior.count=2,plot=TRUE,...)
#	Multidimensional scaling plot of digital gene expression profiles
#	Yunshun Chen, Mark Robinson and Gordon Smyth
#	23 May 2011.  Last modified 11 May 2021.
{
	y <- cpm(x,log=TRUE,prior.count=prior.count)
	return(plotMDS(y,top=top,labels=labels,pch=pch,cex=cex,dim.plot=dim.plot,gene.selection=gene.selection,xlab=xlab,ylab=ylab,plot=plot,...))
}

plotMDS.SummarizedExperiment <- function(x, top=500, labels=NULL, pch=NULL, cex=1, dim.plot=c(1,2), gene.selection="pairwise", xlab=NULL, ylab=NULL, prior.count=2, plot=TRUE, ...)
#	Created 03 April 2020.  Last modified 11 May 2021.
{
	x <- SE2DGEList(x)
	plotMDS.DGEList(x, top=top, labels=labels, pch=pch, cex=cex, dim.plot=dim.plot, gene.selection=gene.selection, xlab=xlab, ylab=ylab, prior.count=prior.count, plot=plot, ...)
}
