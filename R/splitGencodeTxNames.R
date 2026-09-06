splitGencodeTxNames <- function(x, remove.version.numbers = TRUE)
#	Gencode transcript names consist of vertical-bar-separated values,
#	i.e., different annotations pasted together with "|" delimiters.
#	Here we split the separate annotations into the columns of a matrix.
#	Gordon Smyth
#	Created 26 Jun 2026. Last modified 6 Sep 2026.
{
#	Split on 
	A <- strsplit2(x,split="\\|")

#	Remove version numbers from Ensembl transcript and gene IDs
	if(remove.version.numbers) {
		A[,1] <- strsplit2(A[,1],split="\\.")[,1]
		A[,2] <- strsplit2(A[,2],split="\\.")[,1]
	}

	colnames(A) <- c("EnsT","EnsG","HavG","HavT","TxName","GeneName","AnnLength","Type")
	A
}
