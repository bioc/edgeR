\name{PCList}
\alias{PCList}

\title{Create a PCList Object}

\description{
Assembles a \code{PCList} object from its components, especially the paired counts as a couple of matrices.
}

\usage{
PCList(counts, counts2, samples = NULL, group = NULL, genes = NULL, \dots)
}

\arguments{
  \item{counts}{numeric matrix containing sequence read counts, with rows corresponding to genes (genomic features) and columns to samples. Negative values or NAs are not allowed.}
  \item{counts2}{numeric matrix containing sequence read counts, with rows corresponding to genes (genomic features) and columns to samples. Negative values or NAs are not allowed.}
  \item{samples}{data.frame containing sample information, with a row for each sample. This data.frame will be appended to the \code{samples} component of the PCList object.}
  \item{group}{vector or factor giving the experimental group or treatment condition for each sample. Defaults to a single group.}
  \item{genes}{data.frame containing gene annotation.}
  \item{\dots}{other arguments are not currently used.}
}

\details{
Assembles a PCList object from its components.
The only compulsory argument is the paired counts.
}

\value{A \code{\link[edgeR:PCList-class]{PCList}} object.}

\author{edgeR team. Originally created by Lizhong Chen.}

\seealso{\code{\link[edgeR:PCList-class]{PCList-class}}}

\examples{
ngenes <- 100
nsamples <- 4
Counts1 <- matrix(rnbinom(ngenes*nsamples,mu=5,size=10),ngenes,nsamples)
Counts2 <- matrix(rnbinom(ngenes*nsamples,mu=5,size=10),ngenes,nsamples)
rownames(Counts1) <- rownames(Counts2) <- 1:ngenes
colnames(Counts1) <- colnames(Counts2) <- paste0("S",1:4)
Group <- gl(2,2)
Genes <- data.frame(Symbol=paste0("Gene",1:ngenes))
y <- PCList(counts=Counts1, counts2=Counts2, group=Group, genes=Genes)
dim(y)
colnames(y)
y$samples
show(y)
}

\concept{edgeR classes}
