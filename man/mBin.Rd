\name{mbin}
\alias{mbin}
\alias{mbinOneWay}
\alias{mbinIWLS}

\title{Fit Binomial Generalized Linear Models to Multiple Response Vectors: Low Level Functions}

\description{
Fit the same logit-link binomial generalized linear model (GLM) to each row of a matrix of counts.
}

\usage{
mbinOneWay(y, z, design = NULL, group = NULL, offset = NULL, weights = NULL,
        maxit = 50, tol = 1e-8, nthreads = 1L)
mbinIWLS(y, z, design, offset = NULL, weights = NULL,
        maxit = 100, tol = 1e-06, nthreads = 1L)
}

\arguments{
\item{y}{numeric matrix containing the counts.  Rows for genes and columns for libraries.}
\item{z}{numeric matrix containing the paired counts.  Rows for genes and columns for libraries.}

\item{design}{numeric matrix giving the design matrix of the GLM.
Assumed to be full column rank.
This is a required argument for \code{mbinIWLS} and \code{code{designAsFactor}}.
For \code{mBinOneWay}, it defaults to \code{model.matrix(~0+group)} if \code{group} is specified and otherwise to \code{model.matrix( ~1)}.
}

\item{group}{factor giving group membership for oneway layout.
If both \code{design} and \code{group} are specified, then they must agree in terms of \code{designAsFactor}.
If \code{design=NULL}, then a group-means design matrix is implied.}

\item{offset}{numeric vector or matrix giving the offset that is to be included in the log-linear model predictor.  Can be a scalar, a vector of length equal to the number of libraries, or a matrix of the same size as \code{y}.}

\item{weights}{numeric vector or matrix of non-negative quantitative weights.
Can be a vector of length equal to the number of libraries, or a matrix of the same size as \code{y}.}

\item{tol}{numeric scalar giving the convergence tolerance. For \code{mbinOneWay}, convergence is judged successful when the step size falls below \code{tol} in absolute size.}

\item{maxit}{integer giving the maximum number of iterations for the iterated weighted least square algorithm. The iteration will be stopped when this limit is reached even if the convergence criterion hasn't been satisfied.}

\item{nthreads}{integer giving the number of OpenMP threads to use for the genewise fitting. Defaults to 1 (serial).}
}

\details{
These functions are low-level work-horses used by higher-level functions in the edgeR package, especially by \code{binFit}.
Most users will not need to call these functions directly.

The functions \code{mbinOneWay} and \code{mbinIWLS} all fit a binomial GLM to each row of \code{y} and \code{z}.
The row-wise GLMS all have the same design matrix but possibly different offsets and weights.
These functions are all low-level in that they operate on atomic objects (numeric matrices and vectors).

\code{mbinOneWay} fits a oneway layout to each response vector.
It treats the libraries as belonging to a number of groups and fits each group.

\code{mbinIWLS} fits an arbitrary binomial model to each response vector. It implements a standard iterated weighted least square algorithm. The main computation is implemented in C.
}

\value{

\code{mbinOneWay} produces a list with the following components: 
	\item{coefficients}{matrix of estimated coefficients for the models. Rows correspond to rows of \code{y} and columns to columns of \code{design}.}
	\item{fitted.values}{matrix of fitted values. Of same dimensions as \code{y}.}
	\item{deviance}{numeric vector of residual deviances.}

\code{mbinIWLS} produces a list with the following components: 
	\item{coefficients}{matrix of estimated coefficients for the linear models.}
	\item{fitted.values}{matrix of fitted values.}
	\item{deviance}{numeric vector of residual deviances.}
	\item{iter}{number of iterations used.}
	\item{fail}{logical vector indicating genes for which the maximum damping was exceeded before convergence was achieved.}
}

\author{Gordon Smyth, Lizhong Chen.}

\seealso{
Most users will call either \code{\link{binFit}}, the higher-level function offering more object-oriented GLM modeling of Paired Counts data.
}

\concept{Model fit}
