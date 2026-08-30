\name{binFit}
\alias{binFit}
\alias{binFit.PCList}
\alias{binFit.default}

\title{Genewise Binomial Generalized Linear Models}

\description{Fit a binomial generalized log-linear model to the paired counts for each gene.}

\usage{
\method{binFit}{default}(y, z, design=NULL, weights=NULL, offset=NULL, prior.count=1,
                         nthreads=1L, \dots)
\method{binFit}{PCList}(y, design=NULL, \dots)
}

\arguments{
\item{y}{a matrix of counts or a \code{PCList} object. Rows represent genes and columns represent samples.}

\item{z}{a matrix of counts when \code{y} is a matrix. Rows represent genes and columns represent samples.}

\item{design}{numeric matrix giving the design matrix for the genewise linear models.
Must be of full column rank.
Defaults to a single column of ones, equivalent to treating the columns as replicate libraries.}

\item{weights}{positive prior weights for the GLM fits. Can be a single value, or a vector of length \code{ncol(y)}, or a matrix of the same size as \code{y}. If \code{NULL}, will be set to unity for all observations.}

\item{offset}{offsets for the log-linear models.  Can be a single value, or a vector of length \code{ncol(y)}, or a matrix of the same size as \code{y}. If \code{NULL}, it will be set to zero for all observations.}

\item{prior.count}{average prior count to be added to each observation to shrink the estimated log-odds-changes towards zero.}

\item{nthreads}{integer giving the number of OpenMP threads to use for the genewise computations. Defaults to 1 (serial).}

\item{\dots}{other arguments are passed to lower level fitting functions.}
}

\value{
An object of class \code{DGEBIN} (see \code{\link{DGEBIN-class}}) containing the following components:
\item{coefficients}{numeric matrix of estimated coefficients from the glm fits. Coefficients are on the natural log scale. Matrix has dimensions \code{nrow(y)} by \code{ncol(design)}.}
\item{unshrunk.coefficients}{numeric matrix of estimated coefficients from the unshrunk glm fits. Present only when \code{prior.count>0}.}
\item{fitted.values}{matrix of fitted values from glm fits, same number of rows and columns as \code{y}.}
\item{deviance}{numeric vector of deviances, one for each gene.}
\item{df.residual}{numeric vector of residual degrees of freedom, one for each gene.}
\item{counts}{the first matrix of counts as input.}
\item{counts2}{the second matrix of counts as input.}
\item{design}{design matrix as input.}
\item{weights}{matrix of weights as input.}
\item{nthreads}{integer giving the number of OpenMP threads used for the fit, as specified by the \code{nthreads} argument. It is stored on the fitted object so that downstream tests reuse the same number of threads.}
When \code{y} is a \code{PCList}, the components \code{samples}, \code{genes} and \code{AveLogCPM} are also attached.
}

\details{
Implements binomial generalized linear model (GLM) methods developed by Lizhong Chen.
Specifically, \code{binFit} fits genewise binomial GLMs, all with the same design matrix but possibly different offsets and weights.
When the design matrix defines a one-way layout, or can be re-parametrized to a one-way layout, the GLMs are fit very quickly using \code{\link{mbinOneWay}}.
Otherwise the default fitting method, implemented in \code{\link{mbinIWLS}}, uses a standard iterated weighted least square algorithm.
}

\author{Lizhong Chen, Gordon Smyth}

\seealso{
Low-level computations are done by \code{\link{mbinOneWay}} or \code{\link{mbinIWLS}}.

\code{\link{binQLFit}} adds empirical Bayes quasi-dispersion moderation to the binomial fit, and \code{\link{binQLFTest}} performs quasi-likelihood F-tests on the fitted model.
\code{\link{topTags}} displays the results.
}

\examples{
ngenes <- 100
Counts1 <- matrix(rpois(ngenes*4,lambda=20),ngenes,4)
Counts2 <- matrix(rpois(ngenes*4,lambda=20),ngenes,4)
y <- PCList(Counts1, Counts2, group=gl(2,2))
fit <- binFit(y)
}

\concept{Model fit}
\concept{Differential expression}
