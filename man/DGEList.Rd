\name{DGEList}
\alias{DGEList}
\alias{DGEList.default}
\alias{DGEList.data.frame}

\title{
Create a DGEList object
}

\description{
Assembles a \code{DGEList} object from its components, especially the table counts as a matrix or data.frame.
}

\usage{
\method{DGEList}{default}(counts, lib.size = NULL, norm.factors = NULL,
    samples = NULL, group = NULL, genes = NULL,
    remove.zeros = FALSE, \dots)
\method{DGEList}{data.frame}(counts, lib.size = NULL, norm.factors = NULL,
    samples = NULL, group = NULL, genes = NULL,
    remove.zeros = FALSE, annotation.columns = NULL, \dots)
}

\arguments{
  \item{counts}{numeric matrix or data.frame containing sequence read counts, with rows corresponding to genes (genomic features) and columns to samples. Negative values or NAs are not allowed.}
  \item{lib.size}{numeric vector of library sizes (sequencing depths) for the samples. Defaults to \code{colSums(counts)}.}
  \item{norm.factors}{numeric vector of normalization factors that modify the library sizes. Defaults to a vector of ones.}
  \item{samples}{data.frame containing sample information, with a row for each sample. This data.frame will be appended to the \code{samples} component of the DGEList object.}
  \item{group}{vector or factor giving the experimental group or treatment condition for each sample. Defaults to a single group.}
  \item{genes}{data.frame containing gene annotation.}
  \item{remove.zeros}{logical, whether to remove rows that have 0 total count.}
  \item{annotation.columns}{specify columns of \code{counts} that contain gene annotation rather than counts. Can be a vector of column numbers, or a vector of column names, or a logical vector.}
  \item{\dots}{other arguments are not currently used.}
}

\details{
Assembles a DGEList object from its components.
The only compulsory argument is the table of counts.

Normally, \code{counts} is a numeric matrix of counts but a data.frame is also allowed.
If the \code{counts} is a data.frame, then the columns of the data.frame containing gene IDs or other gene annotation can be specified by \code{annotation.columns}, and the other columns are assumed to contain sequence read counts.
If \code{annotation.columns} is not specified, then the function will check for non-numeric columns of \code{counts}
and will attempt to set the leading columns up to the last non-numeric column as annotation.
}

\value{A \code{\link[edgeR:DGEList-class]{DGEList}} object.}

\author{edgeR team. Originally created by Mark Robinson.}

\seealso{\code{\link[edgeR:DGEList-class]{DGEList-class}}}

\examples{
ngenes <- 100
nsamples <- 4
Counts <- matrix(rnbinom(ngenes*nsamples,mu=5,size=10),ngenes,nsamples)
rownames(Counts) <- 1:ngenes
colnames(Counts) <- paste0("S",1:4)
Group <- gl(2,2)
Genes <- data.frame(Symbol=paste0("Gene",1:ngenes))
y <- DGEList(counts=Counts, group=Group, genes=Genes)
dim(y)
colnames(y)
y$samples
show(y)
}

\concept{edgeR classes}
