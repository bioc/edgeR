\name{dim}
\alias{dim.DGEList}
\alias{dim.PCList}
\alias{dim.DGEExact}
\alias{dim.TopTags}
\alias{dim.DGEGLM}
\alias{dim.DGEBIN}
\alias{dim.DGELRT}
\title{Retrieve the Dimensions of a DGEList, PCList, DGEExact, DGEGLM, DGELRT or TopTags Object}
\description{
Retrieve the number of rows (genes) and columns (libraries) for a DGEList, PCList, DGEExact or TopTags Object.
}
\usage{
\method{dim}{DGEList}(x)
}
\arguments{
  \item{x}{an object of class \code{DGEList}, \code{PCList}, \code{DGEExact}, \code{TopTags}, \code{DGEGLM} or \code{DGELRT}.}
}
\details{
Digital gene expression data objects share many analogies with ordinary matrices in which the rows correspond to genes and the columns to arrays.
These methods allow one to extract the size of microarray data objects in the same way that one would do for ordinary matrices.

A consequence is that row and column commands \code{nrow(x)}, \code{ncol(x)} and so on also work.
}
\value{
Numeric vector of length 2.
The first element is the number of rows (genes) and the second is the number of columns (libraries).
}
\author{Gordon Smyth, Davis McCarthy}
\seealso{
  \code{\link[base]{dim}} in the base package.

  \code{\link{DGEList-class}} and \code{\link{PCList-class}} describe the data classes used in edgeR.
}
\examples{
Counts <- matrix(rpois(40,lambda=10),20,2)
y <- DGEList(Counts)
dim(y)
nrow(y)
ncol(y)
}

\concept{edgeR classes}
