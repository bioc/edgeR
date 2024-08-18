\name{normalizeBetweenArrays.DGEList}
\alias{normalizeBetweenArrays.DGEList}
\title{Apply Microarray Normalization to DGEList}

\description{
Applies microrray-style normalization to a DGEList by modifying the offset matrix.
}

\usage{
normalizeBetweenArrays.DGEList(object, method = "cyclicloess", cyclic.method = "affy", \dots)
}

\arguments{
  \item{object}{\code{DGEList} object.}
  \item{method}{character string specifying the normalization method to be used.
  Choices are \code{"none"}, \code{"scale"}, \code{"quantile"} or \code{"cyclicloess"}.}
  \item{cyclic.method}{character string indicating the variant of \code{normalizeCyclicLoess} to be used if \code{method=="cyclicloess"}, see \code{\link{normalizeCyclicLoess}} for possible values.}
  \item{\dots}{other arguments are passed to \code{normalizeQuantiles} or \code{normalizeCyclicLoess}}
}

\details{
This function applies microarray-style normalization methods to a \code{DGEList} object by setting the offset matrix appropriately.
The original counts remain unchanged.

The function \code{cpm()} to compute log2CPM values, then uses \code{limma::normalizeBetweenArrays} to apply the specified normalization method to the log2CPM values, then the difference between the normalized and unnormalized log2CPM matrices is transferred to the offset matrix.
As usual, the \code{offset} entries represent observation-specific normalized loge library sizes.

This function is an alternative to \code{normLibSizes}.
If \code{normLibSizes} has been called previously on the same object, then the results will be ignored.
}

\value{
The input \code{object} is returned with an appropriately normalized \code{offset} component.
}

\author{Gordon Smyth}

\seealso{
\code{\link{normalizeBetweenArrays}}.
}

\examples{
ngenes <- 100
nsamples <- 4
Counts <- matrix(rnbinom(ngenes*nsamples,mu=5,size=10),ngenes,nsamples)
rownames(Counts) <- 1:ngenes
colnames(Counts) <- paste0("S",1:4)
y <- DGEList(counts=Counts)
y <- normalizeBetweenArrays.DGEList(y)
head(y$offset)
}

\keyword{normalization}
