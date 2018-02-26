rowsum.DGEList <- function (x, group, reorder = FALSE, na.rm=FALSE, ...)
#	Sum counts by groups of rows/genes
#	Gordon Smyth
#	22 Feb 2018
{
	isdup <- duplicated(group)
	x2 <- x[!isdup,]
	x2$counts <- rowsum(x$counts,group=group,reorder=FALSE,na.rm=na.rm,...)
	if(reorder) {
		o <- order(row.names(x2))
		x2 <- x2[o,]
	}
	x2
}
