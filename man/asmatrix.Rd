\name{as.matrix}
\alias{as.matrix.DGEList}
\alias{as.matrix.DGEGLM}
\alias{as.matrix.DGEBIN}
\title{Collapse a DGEList, DGEGLM, or DGEBIN Object to a Matrix}
\description{
Coerce a digital gene expression object into a numeric matrix by extracting the count values.
}
\usage{
\method{as.matrix}{DGEList}(x,\dots)
}
\arguments{
  \item{x}{an object of class \code{DGEList}, \code{DGEGLM}, or \code{DGEBIN}.}
  \item{\dots}{additional arguments, not used for these methods.}
}
\details{
This method extracts the matrix of counts from a DGEList or the matrix of coefficients from a DGEGLM or DGEBIN fit.

This involves loss of information, so the original data object is not recoverable.
}
\value{
A numeric matrix.
}
\author{Gordon Smyth}
\seealso{
  \code{\link{as.matrix}} in the base package or \code{\link[limma:as.matrix.EList]{as.matrix}} in the limma package.
}

\concept{Classes}
