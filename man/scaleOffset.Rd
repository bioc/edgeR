\name{scaleOffset}
\alias{scaleOffset}
\alias{scaleOffset.DGEList}
\alias{scaleOffset.default}

\title{Scale Offsets}

\description{Ensures scale of offsets is consistent with library sizes.}

\usage{
\method{scaleOffset}{DGEList}(y, offset, \dots)
\method{scaleOffset}{default}(y, offset, \dots)
}
\arguments{ 
\item{y}{a matrix of counts, or a \code{DGEList} object, or a numeric vector containing library sizes.}
\item{offset}{a numeric vector, matrix, or CompressedMatrix of offsets to be scaled. If a vector, its length must equal the length of \code{y} or the number of columns of \code{y}. If a matrix, its dimension must equal the dimension of \code{y}.}
\item{\dots}{other arguments are not currently used.}
}

\value{
\code{scaleOffset.default} returns a numeric vector if \code{offset} is a vector, a matrix if \code{offset} is a matrix or a \code{CompressedMatrix} object if \code{offset} is a \code{CompressedMatrix}.
\code{scaleOffset.DGEList} computes the scaled offsets and stores them in the \code{offset} component of the input \code{DGEList} object.
}

\details{
This function makes the scale of the offsets is consistent with the library sizes by making mean offset for each gene the same as the mean log-library size.
The length or dimensions of \code{offset} should be consistent with the number of libraries in \code{y}.

The function can be used to ensure that an offset matrix input to edgeR but created by external software is consistent with edgeR's definition of an offset matrix.
It ensures that the each offset entry can be interpreted as the normalized effective log library size for that observation.

This function is not currently called by any other edgeR function, but can be used optionally when inputing an externally created offset matrix to edgeR.
}

\author{Aaron Lun, Yunshun Chen}

\examples{
y <- matrix(rnbinom(40,size=1,mu=100),10,4)
offset <- rnorm(4)
scaleOffset(y, offset)
}

\concept{Normalization}

