plotMDS.PCList <- function (x,top=500,labels=NULL,pch=NULL,cex=1,dim.plot=c(1,2),gene.selection="pairwise",xlab=NULL,ylab=NULL,prior.count=2,plot=TRUE,var.explained=TRUE,...)
#	Multidimensional scaling plot of paired count data
#	Lizhong Chen and Gordon Smyth
#	24 Nov 2025.  Last modified 24 Nov 2025.
{
	# M-values
	M   <- log2(x$counts + prior.count) - log2(x$counts2 + prior.count)

	# compute MDS plot
	mds <- plotMDS(M,top=top,labels=labels,pch=pch,cex=cex,dim.plot=dim.plot,gene.selection=gene.selection,xlab=xlab,ylab=ylab,plot=FALSE,var.explained=var.explained,...)

	# update axislabel
	mds$axislabel <- "Leading logOR dim"

	if (plot)
		plotMDS(mds, labels = labels, pch = pch, cex = cex, xlab = xlab, ylab = ylab, var.explained = var.explained, ...)
	else mds    
}
