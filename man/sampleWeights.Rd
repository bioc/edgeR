\name{sampleWeights}
\alias{sampleWeights}
\title{Empirical Sample Quality Weights}
\description{
Estimate relative quality weights for each sample from the unit deviances and unit residual df output by glmQLFit.
}
\usage{
sampleWeights(unit.deviance.adj, unit.df.adj, s2 = NULL, iter = 10, verbose = FALSE)
}
\arguments{
  \item{unit.deviance.adj}{numeric matrix of (bias corrected) unit deviances. Rows correspond to genes and columns to samples.}
  \item{unit.df.adj}{numeric matrix of (leverage and bias corrected) unit degrees of freedom.}
  \item{s2}{numeric vector of genewise quasi-dispersions. Optional argument, usually left \code{NULL}.}
  \item{iter}{number of iterations.}
  \item{verbose}{logical, if \code{TRUE} then progress information is output for each iteration.}
}
\details{
Estimates an empirical quality weight for each sample from the matrix of adjusted unit deviances and adjusted unit df
from \code{glmQLFit} with \code{legacy=FALSE} and \code{keep.unit.mat=TRUE}.
The quasi-dispersion of each observation is assumed to be equal to \code{s2_g / w_i}.

The vector of sample weights can be entered as the \code{weights} argument for subsequent \code{glmQLFit} fits.
}
\value{
Numeric vector of length \code{ncol(fit)} containing the sample weights.
}
\author{Gordon Smyth}
\seealso{
This function uses output from \code{\link{glmQLFit}} and is analogous to \code{\link[limma]{arrayWeights}} in the limma package.
The \code{\link[limma]{voomWithQualityWeights}} and \code{\link{voomLmFit}} functions also estimate sample weights for count data.
}
\examples{
# Sample 2 has higher dispersion and lower weight
y <- matrix(rpois(1000*5,lambda=8),1000,5)
y[,2] <- 2 * y[,2]

# Usually the NB dispersion is estimated by glmQLFit(),
# but here we set it to the true value for Poisson counts to illustrate the weights in pure form
fit <- glmQLFit(y, dispersion=0, keep.unit.mat=TRUE)
w <- sampleWeights(fit$unit.deviance.adj, fit$unit.df.adj)
barplot(w, main="Sample weights")
}
