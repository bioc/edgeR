splitIntoGroupsPseudo <- function(pseudo,group,pair) 
#	Extract the data for specified two groups from a matrix of pseudocounts pair <- levels(as.factor(pair))
#	Written by Davis McCarthy, Feb 2009, idea suggested by Mark Robinson
{
	y1 <- pseudo[, group==pair[1], drop=FALSE]
	y2 <- pseudo[, group==pair[2], drop=FALSE]
	list(y1=y1, y2=y2)
}
