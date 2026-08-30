\title{Multidimensional Scaling Plot of Distances between M-values from the Paired Count Data}
\name{plotMDS.PCList}
\alias{plotMDS.PCList}
\description{
Plot samples on a two-dimensional scatterplot so that distances on the plot approximate the M-values differences between the samples.
}
\usage{
\method{plotMDS}{PCList}(x, top = 500, labels = NULL, pch = NULL, cex = 1,
       dim.plot = c(1,2), gene.selection = "pairwise", xlab = NULL, ylab = NULL,
       prior.count = 2, plot = TRUE, var.explained = TRUE, \dots)
}
\arguments{
  \item{x}{a \code{PCList} object.}
  \item{top}{number of top features used to calculate pairwise distances.}
  \item{labels}{character vector of sample names or labels. If \code{x} has no column names, then it defaults to the index of the samples.}
  \item{pch}{plotting symbol or symbols. See \code{\link{points}} for possible values. Ignored if \code{labels} is non-\code{NULL}.}
  \item{cex}{numeric vector of plot symbol expansions. See \code{\link[graphics]{text}} for possible values.}
  \item{dim.plot}{which two dimensions should be plotted, numeric vector of length two.}
  \item{gene.selection}{character, \code{"pairwise"} to choose the top genes separately for each pairwise comparison between the samples, or \code{"common"} to select the same genes for all comparisons.}
  \item{xlab}{x-axis label.}
  \item{ylab}{y-axis label.}
  \item{prior.count}{average prior count added to each count before forming the M-values (log2 odds), to shrink the values for low counts towards zero.}
  \item{plot}{logical. If \code{TRUE} then a plot is created on the current graphics device.}
  \item{var.explained}{logical. If \code{TRUE} then the percentage variation explained is included in the axis labels.}
  \item{\dots}{any other arguments are passed to \code{plot}.}
}

\details{
The paired count data is converted to M-values (log2 odds) and passed to the limma \code{plotMDS} function.
This method calculates distances between samples based on log2 odds ratios.
See the \code{\link[limma:plotMDS]{plotMDS help page}} for details.

The number of features (\code{top}) chosen for this exercise should roughly correspond to the number of features with large log2 odds ratios.
The default setting of 500 features is widely effective and suitable for routine use, but a smaller value might be chosen when the samples are distinguished by a specific focused molecular pathway.
Very large values (greater than 1000) are not usually so effective.
}

\value{
An object of class \code{\link[limma:plotMDS]{MDS}} is invisibly returned and (if \code{plot=TRUE}) a plot is created on the current graphics device.
}

\author{Lizhong Chen and Gordon Smyth}

\seealso{
\code{\link[limma]{plotMDS}}, \code{\link{cmdscale}}, \code{\link{as.dist}}
}

\examples{
ngenes <- 100
Counts1 <- matrix(rpois(ngenes*4,lambda=10),ngenes,4)
Counts2 <- matrix(rpois(ngenes*4,lambda=10),ngenes,4)
y <- PCList(Counts1, Counts2, group=gl(2,2))
plotMDS(y)
}

\concept{Data exploration}

