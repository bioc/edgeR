\name{binQLFTest}
\alias{binQLFTest}

\title{Quasi-Likelihood F-Tests for Binomial Models}

\description{Conduct genewise quasi-likelihood F-tests for a given coefficient or coefficient contrast, for binomial models fitted by \code{\link{binQLFit}}.}

\usage{
binQLFTest(glmfit, coef = ncol(glmfit$design), contrast = NULL, or = 1, lor = NULL,
        upshot = TRUE)
}

\arguments{
\item{glmfit}{a \code{DGEBIN} object, usually output from \code{\link{binQLFit}}. Must contain the empirical Bayes quasi-dispersions (\code{s2.post}).}

\item{coef}{integer or character index vector indicating which coefficients of the binomial linear model are to be tested equal to zero. Ignored if \code{contrast} is not \code{NULL}.}

\item{contrast}{numeric vector or matrix specifying one or more contrasts of the linear model coefficients to be tested equal to zero.}

\item{or}{scalar or numeric vector of minimum odds ratio below which changes are not considered scientifically meaningful.}

\item{lor}{scalar or numeric vector of minimum absolute log2 odds ratio below which changes are not considered scientifically meaningful. Defaults to \code{log2(or)}. If specified then takes precedence over \code{or}.}

\item{upshot}{logical. If \code{TRUE}, then a more powerful test is conducted assuming that the true log2 odds ratios are uniformly distributed over the null interval from \code{-lor} to \code{lor}. If \code{FALSE}, then the standard TREAT p-values are computed assuming the worst-case scenario of true log2 odds ratio equal to \code{lor}.}
}

\details{
\code{binQLFTest} is the binomial analogue of \code{\link{glmQLFTest}} and is typically used after \code{\link{binQLFit}}.
It conducts empirical Bayes quasi-likelihood F-tests for the binomial generalized linear models fitted genewise to paired counts (for example methylated and unmethylated read counts) by \code{\link{binQLFit}}.

The recommended workflow is to build a \code{\link{PCList}} object, fit the models with \code{\link{binQLFit}}, then test with \code{binQLFTest}:
\code{binQLFit} attaches the average log2-CPM (\code{AveLogCPM}) and the moderated quasi-dispersions that \code{binQLFTest} requires.

The reported \code{logOR} is the log2 odds ratio for the tested coefficient or contrast, equal to the estimated coefficient divided by \code{log(2)}.
When \code{binQLFit} is run with \code{prior.count>0} (the default), \code{logOR} is reported on the shrunk coefficient scale, whereas the F-test and its deviances are computed from the unshrunk maximum-likelihood fit.

When a nonzero threshold is specified via \code{or} or \code{lor}, a TREAT-style test relative to the threshold is performed.
With \code{upshot=TRUE} the more powerful UPSHOT p-values are computed by Gauss quadrature over the null interval.
}

\value{
An object of class \code{DGELRT} with the same components as \code{glmfit} plus the following:
\item{table}{data frame with the same rows as \code{glmfit} and the following columns:
  \code{logOR}, the log2 odds ratio for the tested coefficient or contrast;
  \code{logCPM}, the average log2-counts per million;
  \code{F}, the quasi-likelihood F-statistic;
  \code{PValue}, the associated p-value.}
\item{comparison}{character string describing the coefficient or contrast being tested.}
\item{df.test}{numerator degrees of freedom for the F-test, equal to the number of coefficients tested.}
\item{df.total}{numeric vector of denominator (total) degrees of freedom for the F-test.}
}

\author{Lizhong Chen, Gordon Smyth}

\seealso{
\code{\link{binQLFit}} fits the binomial models that are provided to \code{binQLFTest}, and \code{\link{PCList}} constructs the paired-count data object.

\code{\link{topTags}} displays the results in a convenient form.

\code{\link{glmQLFTest}} performs the analogous quasi-likelihood F-tests for negative binomial models.
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
