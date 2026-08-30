\name{dimnames}
\alias{dimnames.DGEList}
\alias{dimnames.PCList}
\alias{dimnames.DGEExact}
\alias{dimnames.DGEGLM}
\alias{dimnames.DGEBIN}
\alias{dimnames.DGELRT}
\alias{dimnames.TopTags}
\alias{dimnames<-.DGEList}
\alias{dimnames<-.PCList}
\alias{dimnames<-.DGEExact}
\alias{dimnames<-.DGEGLM}
\alias{dimnames<-.DGEBIN}
\alias{dimnames<-.DGELRT}

\title{Retrieve the Dimension Names of a DGE Object}
\description{
Retrieve the dimension names of a digital gene expression data object.
}
\usage{
\method{dimnames}{DGEList}(x)
\method{dimnames}{DGEList}(x) <- value
}
\arguments{
  \item{x}{an object of class \code{DGEList}, \code{PCList}, \code{DGEExact}, \code{DGEGLM}, \code{DGELRT} or \code{TopTags}.}
  \item{value}{a possible value for \code{dimnames(x)}, see \code{\link{dimnames}}.}
}
\details{
The dimension names of a DGE data object are the same as those of the most important component of that object.

Setting dimension names is supported for \code{DGEList}, \code{PCList}, \code{DGEExact}, \code{DGEGLM} and \code{DGELRT} objects.

A consequence of these methods is that \code{rownames}, \code{colnames}, \code{rownames<-} and \code{colnames<-} will also work as expected on any of the above object classes.
}
\value{
Either \code{NULL} or a list of length 2.
If a list, its components are either \code{NULL} or a character vector the length of the appropriate dimension of \code{x}.
}
\author{Gordon Smyth}
\seealso{
  \code{\link[base]{dimnames}} in the base package.
}

\examples{
Counts1 <- matrix(rpois(20,lambda=10),10,2)
Counts2 <- matrix(rpois(20,lambda=10),10,2)
y <- PCList(Counts1, Counts2)
rownames(y) <- paste0("Locus",1:10)
colnames(y) <- c("A","B")
dimnames(y)
}

\concept{edgeR classes}