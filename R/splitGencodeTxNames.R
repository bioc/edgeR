splitGencodeTxNames <- function(x, remove.version=TRUE)
#	Split GENCODE transcript IDs into an 8-column character matrix.
#	Gencode transcript names consist of vertical-bar-separated values,
#	i.e., different annotations pasted together with "|" delimiters.
#	Gordon Smyth
#	Created 26 Jun 2026. Last modified 22 Oct 2026.
{
#	Split pipe-delimited IDs. True Gencode tx IDs should yield 8 columns.
	x <- as.character(x)
	A <- strsplit2(x, split="|", fixed=TRUE)

#	Remove version numbers from Ensembl transcript and gene IDs
	j <- min(2L,ncol(A))
	if(remove.version) {
		A[,1:j] <- sub("\\.[0-9]+", "", A[,1:j])
	}

#	Abbreviated column names: Ens for Ensembl, Hav for Havana, Tx for Transcript, Type for Biotype
	GencodeFields <- c("EnsT","EnsG","HavG","HavT","TxName","GeneName","AnnLength","Type")
	if(ncol(A)==8L) {
		colnames(A) <- GencodeFields
	} else {
		if(ncol(A) < 8L) {
			colnames(A) <- GencodeFields[1:ncol(A)]
		} else {
			colnames(A) <- c(GencodeFields,paste0("V",9:ncol(A)))
		}
	}
	A
}

collapseGencodeTx2Gene <- function(txann)
#	Collapse output from splitGenecodeTxNames() to gene-level data.fame.
#	Gordon Smyth
#	Created 8 Oct 2026. Last modified 9 Oct 2026.
{
#	Unique gene IDs
	EnsG <- txann[,"EnsG"]
	d <- duplicated(EnsG)
	EnsGu <- EnsG[!d]

#	Abbreviate snake_case tx types
	TxType <- txann[,"Type"]
	HasUnderscore <- grepl("_",TxType,fixed=TRUE)
	TxType[HasUnderscore] <- gsub("([^_]{3})[^_]*","\\1",TxType[HasUnderscore])

#	Collapse unique tx types into a composite type for each gene
	GeneType <- vapply(
		split(TxType, factor(EnsG,levels=EnsGu)),
		function(x) paste(sort(unique(x)), collapse="|"),
	 	character(1)
	 )

#	Further abbreviation
	GeneType <- gsub("pro_cod_CDS_not_def","CDS_not_def",GeneType)
	GeneType <- gsub("pro_cod_LoF","LoF",GeneType)

#	Assemble output data.frame
	out <- as.data.frame(txann[!d,"GeneName",drop=FALSE])
	row.names(out) <- EnsGu
	out$TxTypes <- GeneType
	out$NPCTx <- drop(rowsum(as.numeric(txann[,"Type"]=="protein_coding"),EnsG,reorder=FALSE))
	out
}


