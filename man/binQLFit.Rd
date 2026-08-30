\name{binQLFit}
\alias{binQLFit}
\alias{binQLFit.PCList}
\alias{binQLFit.default}

\title{Quasi-Likelihood Binomial Generalized Linear Models}

\description{Fit a binomial generalized log-linear model to the paired counts for each gene and estimate the empirical Bayes moderated quasi-dispersions.}

\usage{
\method{binQLFit}{default}(y, z, design=NULL, weights=NULL, offset=NULL, AveLogCPM=NULL,
                           abundance.trend=TRUE, covariate.trend=NULL, robust=TRUE,
                           prior.count=1, nthreads=1L, \dots)
\method{binQLFit}{PCList}(y, design=NULL, abundance.trend=TRUE, covariate.trend=NULL, robust=TRUE, \dots)
}

\arguments{
\item{y}{a matrix of counts or a \code{PCList} object. Rows represent genes and columns represent samples.}

\item{z}{a matrix of counts when \code{y} is a matrix. Rows represent genes and columns represent samples.}

\item{design}{numeric matrix giving the design matrix for the genewise linear models.
Must be of full column rank.
Defaults to a single column of ones, equivalent to treating the columns as replicate libraries.}

\item{weights}{positive prior weights for the GLM fits. Can be a single value, or a vector of length \code{ncol(y)}, or a matrix of the same size as \code{y}. If \code{NULL}, will be set to unity for all observations.}

\item{offset}{offsets for the log-linear models.  Can be a single value, or a vector of length \code{ncol(y)}, or a matrix of the same size as \code{y}. If \code{NULL}, it will be set to zero for all observations.}

\item{AveLogCPM}{numeric vector giving average log2-counts per million for each row of \code{y}. If \code{NULL}, then will be computed by \code{aveLogCPM(y)}.}

\item{abundance.trend}{logical, whether to allow an abundance trend for the quasi-dispersion prior.}

\item{covariate.trend}{numeric vector of length \code{nrow(y)}. If non-NULL, then will be used instead of \code{AveLogCPM} as the covariate for the trended quasi-dispersion prior.}

\item{robust}{logical, whether to estimate the prior distribution for the quasi-dispersions robustly.}

\item{prior.count}{average prior count to be added to each observation to shrink the estimated log-odds-changes towards zero.}

\item{nthreads}{integer giving the number of OpenMP threads to use for the genewise computations. Defaults to 1 (serial).}

\item{\dots}{other arguments are passed to lower level fitting functions.}
}

\value{
An object of class \code{DGEBIN} (see \code{\link{DGEBIN-class}}) with the raw-fit components produced by \code{\link{binFit}} plus the following:
\item{df.residual.adj}{numeric vector of adjusted residual degrees of freedom for each gene.}
\item{deviance.adj}{numeric vector of adjusted deviances, one for each gene.}
\item{df.prior}{numeric vector of prior degrees of freedom for the quasi-dispersions. Has length \code{nrow(y)} if \code{robust=TRUE}, otherwise length 1.}
\item{s2.prior}{numeric vector giving prior value for the quasi-dispersions. Has length \code{nrow(y)} if \code{abundance.trend=TRUE}, otherwise length 1.}
\item{s2.post}{numeric vector of posterior genewise quasi-dispersions.}
\item{average.ql.dispersion}{average quasi-dispersion, the overall overdispersion for binomial models.}
\item{AveLogCPM}{numeric vector giving average log2-counts per million for each gene.}
}

\details{
\code{binQLFit} is the binomial analogue of \code{\link{glmQLFit}}.
It first fits genewise binomial GLMs (as done by \code{\link{binFit}}), all with the same design matrix but possibly different offsets and weights, then estimates the genewise quasi-dispersions and moderates them towards a common prior by empirical Bayes using \code{\link[limma]{squeezeVar}}.
When the design matrix defines a one-way layout, or can be re-parametrized to a one-way layout, the GLMs are fit very quickly using \code{\link{mbinOneWay}}.
Otherwise the default fitting method, implemented in \code{\link{mbinIWLS}}, uses a standard iterated weighted least square algorithm.

The moderated quasi-dispersions (\code{s2.post}) and adjusted degrees of freedom are required by \code{\link{binQLFTest}} to conduct the quasi-likelihood F-tests.
}

\author{Lizhong Chen, Gordon Smyth}

\seealso{
\code{\link{binFit}} performs the underlying binomial GLM fit without the quasi-dispersion moderation.

\code{\link{binQLFTest}} performs quasi-likelihood F-tests on the fitted model and \code{\link{topTags}} displays the results.

\code{\link{glmQLFit}} performs the analogous quasi-likelihood fit for negative binomial models.
}

\examples{
ngenes <- 100
Counts1 <- matrix(rpois(ngenes*4, lambda=20), ngenes, 4)
Counts2 <- matrix(rpois(ngenes*4, lambda=20), ngenes, 4)
y <- PCList(Counts1, Counts2, group=gl(2,2))
fit <- binQLFit(y)
res <- binQLFTest(fit, coef=2)
topTags(res)
}

\concept{Model fit}
\concept{Differential expression}
