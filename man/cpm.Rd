\name{cpm}
\alias{cpm}
\alias{cpm.DGEList}
\alias{cpm.default}
\alias{rpkm}
\alias{rpkm.DGEList}
\alias{rpkm.default}
\alias{cpmByGroup}
\alias{cpmByGroup.DGEList}
\alias{cpmByGroup.default}
\alias{rpkmByGroup}
\alias{rpkmByGroup.DGEList}
\alias{rpkmByGroup.default}

\title{Counts per Million or Reads per Kilobase per Million}

\description{Compute counts per million (CPM) or reads per kilobase per million (RPKM).}

\usage{
\method{cpm}{DGEList}(y, normalized.lib.sizes = TRUE,
       log = FALSE, prior.count = 2, \dots)
\method{cpm}{default}(y, lib.size = NULL,
       log = FALSE, prior.count = 2, \dots)
\method{rpkm}{DGEList}(y, gene.length = NULL, normalized.lib.sizes = TRUE,
       log = FALSE, prior.count = 2, \dots)
\method{rpkm}{default}(y, gene.length, lib.size = NULL,
       log = FALSE, prior.count = 2, \dots)
\method{cpmByGroup}{DGEList}(y, group = NULL, dispersion = NULL, \ldots)
\method{cpmByGroup}{default}(y, group = NULL,
       dispersion = 0.05, offset = NULL, weights = NULL, \ldots)
\method{rpkmByGroup}{DGEList}(y, group = NULL, gene.length = NULL, dispersion = NULL, \ldots)
\method{rpkmByGroup}{default}(y, group = NULL, gene.length,
       dispersion = 0.05, offset = NULL, weights = NULL, \ldots)
}

\arguments{ 
\item{y}{matrix of counts or a \code{DGEList} object}
\item{normalized.lib.sizes}{logical, use normalized library sizes?}
\item{lib.size}{library size, defaults to \code{colSums(y)}.}
\item{log}{logical, if \code{TRUE} then \code{log2} values are returned.}
\item{prior.count}{average count to be added to each observation to avoid taking log of zero. Used only if \code{log=TRUE}.}
\item{gene.length}{vector of length \code{nrow(y)} giving gene length in bases, or the name of the column \code{y$genes} containing the gene lengths.}
\item{group}{factor giving group membership for columns of \code{y}. Defaults to \code{y$sample$group} for the \code{DGEList} method and to a single level factor for the default method.}
\item{dispersion}{numeric vector of negative binomial dispersions.}
\item{offset}{numeric matrix of same size as \code{y} giving offsets for the log-linear models.  Can be a scalar or a vector of length \code{ncol(y)}, in which case it is expanded out to a matrix.}
\item{weights}{numeric vector or matrix of non-negative quantitative weights.
Can be a vector of length equal to the number of libraries, or a matrix of the same size as \code{y}.}
\item{\dots}{other arguments are not used.}
}

\value{A numeric matrix of CPM or RPKM values.
\code{cpm} and \code{rpkm} produce matrices of the same size as \code{y}.
\code{cpmByGroup} and \code{rpkmByGroup} produce matrices with a column for each level of \code{group}.
If \code{log = TRUE}, then the values are on the log2 scale.}

\details{
CPM or RPKM values are useful descriptive measures for the expression level of a gene.
By default, the normalized library sizes are used in the computation for \code{DGEList} objects but simple column sums for matrices.

If log-values are computed, then a small count, given by \code{prior.count} but scaled to be proportional to the library size, is added to \code{y} to avoid taking the log of zero.

The \code{rpkm} method for \code{DGEList} objects will try to find the gene lengths in a column of \code{y$genes} called \code{Length} or \code{length}.
Failing that, it will look for any column name containing \code{"length"} in any capitalization.

\code{cpmByGroup} and \code{rpkmByGroup} compute group average values on the unlogged scale.
}

\note{
\code{aveLogCPM(y)}, \code{rowMeans(cpm(y,log=TRUE))} and \code{log2(rowMeans(cpm(y))} all give slightly different results.
}

\author{Davis McCarthy, Gordon Smyth}

\seealso{
\code{\link{aveLogCPM}}
}

\examples{
y <- matrix(rnbinom(20,size=1,mu=10),5,4)
cpm(y)

d <- DGEList(counts=y, lib.size=1001:1004)
cpm(d)
cpm(d,log=TRUE)

d$genes <- data.frame(Length=c(1000,2000,500,1500,3000))
rpkm(d)
}
