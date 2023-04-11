\name{Seurat2PB}
\alias{Seurat2PB}

\title{Seurat or SeuratObject class object to pseudo-bulk DGEList}

\description{Given a \code{Seurat} or \code{SeuratObject} data object, create pseudo-bulk samples using the sample and cluster information and return a \code{DGEList} object.}

\usage{
Seurat2PB(object, sample, cluster="seurat_clusters")
}

\arguments{
\item{object}{a \code{Seurat} or \code{SeuratObject} class data object. Must have \code{counts} in its \code{RNA} assay.}
\item{sample}{character specifying the column of \code{object@meta.data} that contains sample information of the single cell data. Must be one of the column names of \code{object@meta.data}.}
\item{cluster}{character specifying the column of \code{object@meta.data} that contains single cell cluster information. Must be one of the column names of \code{object@meta.data}. Default to \code{seurat_clusters}.}
}

\details{
This function takes a \code{Seurat} or \code{SeuratObject} data object as input.
It is assumed that the input data object contains raw RNA-seq read counts of multiple samples as well as cell clustering information (usually from a single cell integration analysis).
The \code{counts} of the \code{RNA} assay of the input data object is first extracted.
A pseudo-bulk count matrix is formed by aggregating the read counts of all the cells within the same \code{cluster} for each \code{sample}.
The pseudo-bulk count matrix is then used as the \code{counts} component of the output \code{DGEList} object.
The \code{cluster} and \code{sample} information of the pseudo-bulk counts is stored in the \code{samples} component.
The row names of the input data object, together with other extra feature information (if any) found in the input, are stored in the \code{genes} component.
}

\value{
A \code{DGEList} object.
}

\author{Yunshun Chen}

\examples{
\dontrun{
ngenes <- 1e3
ncells <- 1e3
sp <- paste0("sample", sample(1:2, ncells, replace=TRUE))
clst <- sample(1:3, ncells, replace=TRUE)

counts <- matrix(rnbinom(ngenes*ncells, mu=10, size=2), ngenes, ncells)
colnames(counts) <- paste0("Cell",1:ncells)
rownames(counts) <- paste0("Gene",1:ngenes)
so <- CreateSeuratObject(counts = counts)
so@meta.data <- cbind(so@meta.data, sample=sp, cluster=clst)

y <- Seurat2PB(so, sample="sample", cluster="cluster")
}}

\concept{edgeR classes}

