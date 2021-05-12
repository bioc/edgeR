\title{Multidimensional scaling plot of distances between digital gene expression profiles}
\name{plotMDS.DGEList}
\alias{plotMDS.DGEList}
\alias{plotMDS.SummarizedExperiment}
\description{
Plot samples on a two-dimensional scatterplot so that distances on the plot approximate the expression differences between the samples.
}
\usage{
\method{plotMDS}{DGEList}(x, top = 500, labels = NULL, pch = NULL, cex = 1,
        dim.plot = c(1,2), gene.selection = "pairwise", xlab = NULL, 
        ylab = NULL, prior.count = 2, plot = TRUE, \dots)
\method{plotMDS}{SummarizedExperiment}(x, top = 500, labels = NULL, pch = NULL, cex = 1,
        dim.plot = c(1,2), gene.selection = "pairwise", xlab = NULL, 
        ylab = NULL, prior.count = 2, plot = TRUE, \dots)
}
\arguments{
  \item{x}{a \code{DGEList} or \code{SummarizedExperiment} object.}
  \item{top}{number of top genes used to calculate pairwise distances.}
  \item{labels}{character vector of sample names or labels. If \code{x} has no column names, then defaults the index of the samples.}
  \item{pch}{plotting symbol or symbols. See \code{\link{points}} for possible values. Ignored if \code{labels} is non-\code{NULL}.}
  \item{cex}{numeric vector of plot symbol expansions. See \code{\link[graphics]{text}} for possible values.}
  \item{dim.plot}{which two dimensions should be plotted, numeric vector of length two.}
  \item{gene.selection}{character, \code{"pairwise"} to choose the top genes separately for each pairwise comparison between the samples, or \code{"common"} to select the same genes for all comparisons.}
  \item{xlab}{x-axis label}
  \item{ylab}{y-axis label}
  \item{prior.count}{average prior count to be added to observation to shrink the estimated log-fold-changes towards zero.}
  \item{plot}{logical. If \code{TRUE} then a plot is created on the current graphics device.}
  \item{\dots}{any other arguments are passed to \code{plot}.}
}

\details{
The function converts the counts to log-counts-per-million using \code{cpm} and to pass these to the limma \code{plotMDS} function.
It calculates distances between samples based on log2 fold changes.
See the \code{\link[limma:plotMDS]{plotMDS help page}} for details.

The number of genes (\code{top}) chosen for this exercise should roughly correspond to the number of differentially expressed genes with materially large fold-changes.
The default setting of 500 genes is widely effective and suitable for routine use, but a smaller value might be chosen for when the samples are distinguished by a specific focused molecular pathway.
Very large values (greater than 1000) are not usually so effective.
}

\value{
An object of class \code{\link[limma:plotMDS]{MDS}} is invisibly returned and (if \code{plot=TRUE}) a plot is created on the current graphics device.
}

\author{Yunshun Chen, Mark Robinson and Gordon Smyth}

\seealso{
\code{\link[limma]{plotMDS}}
}

\examples{
# Simulate DGE data for 1000 genes and 6 samples.
# Samples are in two groups
# First 200 genes are differentially expressed in second group

ngenes <- 1000
nlib <- 6
counts <- matrix(rnbinom(ngenes*nlib, size=1/10, mu=20),ngenes,nlib)
rownames(counts) <- paste("gene",1:ngenes, sep=".")
group <- gl(2,3,labels=c("Grp1","Grp2"))
counts[1:200,group=="Grp2"] <- counts[1:200,group=="Grp2"] + 10
y <- DGEList(counts,group=group)
y <- calcNormFactors(y)

# without labels, indexes of samples are plotted.
col <- as.numeric(group)
mds <- plotMDS(y, top=200, col=col)

# or labels can be provided, here group indicators:
plotMDS(mds, col=col, labels=group)
}

\concept{Data exploration}
