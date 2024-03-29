\name{as.matrix}
\alias{as.matrix.DGEList}
\alias{as.matrix.DGEGLM}
\title{Collapse a DGEList or DGEGLM Object to a Matrix}
\description{
Coerce a digital gene expression object into a numeric matrix by extracting the count values.
}
\usage{
\method{as.matrix}{DGEList}(x,\dots)
}
\arguments{
  \item{x}{an object of class \code{DGEList} or \code{DGEGLM}.}
  \item{\dots}{additional arguments, not used for these methods.}
}
\details{
This method extracts the matrix of counts from a DGEList or the matrix of coefficients from a DGEGLM fit.

This involves loss of information, so the original data object is not recoverable.
}
\value{
A numeric matrix.
}
\author{Gordon Smyth}
\seealso{
  \code{\link{as.matrix}} in the base package or \code{\link[limma:asmatrix]{as.matrix}} in the limma package.
}

\concept{edgeR classes}
