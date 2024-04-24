Seurat2PB <- function(object, sample, cluster="seurat_clusters")
#	Given a 'Seurat' or 'SeuratObject' data object, create pseudo-bulk samples using
#   the sample and cluster information and return a DGEList object
#	Yunshun Chen
#	04 Jan 2023. Last modified 24 Apr 2024.
{
	if(!is(object,"SeuratObject") & !is(object,"Seurat"))
		stop("object must belong to either the SeuratObject or Seurat classes")

	if(!requireNamespace("SeuratObject",quietly=TRUE))
		stop("SeuratObject package required but is not installed (or can't be loaded)")

#	Check 'assays'
    if(packageVersion("SeuratObject") < "5.0.0"){
        counts <- SeuratObject::GetAssayData(object, assay="RNA", slot="counts")
    } else {
        counts <- SeuratObject::GetAssayData(object, assay="RNA", layer="counts")
    }
	if( is.null(counts) ) stop("object doesn't contain raw RNA counts")

#   Check 'meta.data'
    meta <- object@meta.data
    if(! sample %in% names(meta) ) stop("sample information can not be found in meta.data")
    if(! cluster %in% names(meta) ) stop("cluster information can not be found in meta.data")

    sp <- meta[, sample]
    clst <- meta[, cluster]
    if(length(table(sp)) == 1) warning("Only 1 sample found in meta.data. Please check whether sample information is specified correctly.")
    if(length(table(clst)) == 1) warning("Only 1 cluster found in meta.data. Please check whether cluster information is specified correctly.")

#	Check gene information
    genes <- data.frame( gene=rownames(object[["RNA"]]) )
    genes <- cbind( genes, object[["RNA"]][[]])
    
#   Pseudo-bulk counts
    sp_clst <- factor(paste(sp, clst, sep="_cluster"))
    group_mat <- Matrix::sparse.model.matrix(~ 0 + sp_clst)
    colnames(group_mat) <- gsub("^sp_clst", "", colnames(group_mat))
    counts.pb <- counts %*% group_mat

#   Pseudo-bulk sample information
    levels(sp_clst)
    sp.pb <- gsub("_cluster.*$", "", levels(sp_clst))
    clst.pb <- gsub("^.*_cluster", "", levels(sp_clst))
    sample.pb <- data.frame(sample=sp.pb, cluster=clst.pb)

#   DGEList
	DGEList(counts=as.matrix(counts.pb), samples=sample.pb, genes=genes)
}

