\name{loessByCol}
\alias{loessByCol}
\alias{locfitByCol}

\title{Locally Weighted Mean by Column}

\description{Smooth columns of a matrix by non-robust loess curves of degree 0.}

\usage{
loessByCol(y, x=NULL, span=0.5)
locfitByCol(y, x=NULL, weights=1, span=0.5, degree=0, nthreads=1L)
}

\arguments{
\item{y}{numeric matrix of response variables.}

\item{x}{numeric covariate vector of length \code{nrow(y)}, defaults to equally spaced.}

\item{span}{width of the smoothing window, in terms of proportion of the data set. Larger values produce smoother curves.}

\item{weights}{relative weights of each observation, one for each covariate value.}

\item{degree}{degree of local polynomial fit.}

\item{nthreads}{integer giving the number of OpenMP threads to use for fitting the columns. Defaults to 1 (serial). Results do not depend on it.}
}

\value{A list containing a numeric matrix with smoothed columns and a vector of leverages for each covariate value.

\code{locfitByCol} returns a numeric matrix.
}

\details{
Fits a loess curve with degree 0 to each column of the response matrix, using the same covariate vector for each column.
The smoothed column values are tricube-weighted means of the original values.

\code{locfitByCol} returns, for each column, the fitted values of \code{locfit(y[,j] ~ x, weights=weights, alpha=span, deg=degree)} with the other \code{locfit} settings at their defaults.
The fit is computed by C code implementing the \code{locfit} algorithm (Loader, 1999): tricube local regression with nearest-neighbour bandwidths, fitted at the vertices of an adaptive tree and interpolated between them.
The results agree with \code{locfit} up to rounding error.
}

\author{Aaron Lun for \code{loessByCol}, replacing earlier R code by Davis McCarthy. Gordon Smyth for \code{locfitByCol}. Lizhong Chen for new C versions.}

\seealso{
\code{\link{loess}}
}

\examples{
y <- matrix(rnorm(100*3), nrow=100, ncol=3)
head(y)
out <- loessByCol(y)
head(out$fitted.values)
}
