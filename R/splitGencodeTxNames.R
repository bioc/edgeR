splitGencodeTxNames <- function(x)
# Genecode transcript names consist of different annotations pasted together with "|" delimiter.
# Here we split the separate annotations into columns of a matrix.
# Gordon Smyth
# 26 Jun 2026.
{
	A <- strsplit2(x,split="\\|")
#	Remove version number from Ensembl transcript ID
	A[,1] <- strsplit2(A[,1],split="\\.")[,1]
#	Remove version number from Ensembl gene ID
	A[,2] <- strsplit2(A[,2],split="\\.")[,1]
	colnames(A) <- c("EnsT","EnsG","HavG","HavT","TxName","GeneName","AnnLength","Type")
	A
}
