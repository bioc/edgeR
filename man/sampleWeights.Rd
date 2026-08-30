\name{sampleWeights}
\alias{sampleWeights}
\title{Empirical Sample Quality Weights}
\description{
Estimate relative quality weights for each sample from the unit deviances and unit residual df output by glmQLFit.
}
\usage{
sampleWeights(glmfit)
}
\arguments{
  \item{glmfit}{a \code{DGEGLM} object (from \code{glmQLFit}) or a \code{DGEBIN} object (from \code{binQLFit}). The number of OpenMP threads is taken from \code{glmfit$nthreads}.}
}
\details{
Estimates an empirical quality weight for each sample from the matrix of adjusted unit deviances and adjusted unit df
from \code{glmQLFit} (a \code{DGEGLM} object) or \code{binQLFit} (a \code{DGEBIN} object), both with \code{legacy=FALSE}.
The quasi-dispersion of each observation is assumed to be equal to \code{s2_g / w_i}.

The returned sample weights can be entered directly as the \code{weights} argument for subsequent \code{glmQLFit} or \code{binQLFit} fits.
}
\value{
A \code{CompressedMatrix} with the same dimensions as \code{glmfit$counts} holding the per-sample weights, repeated across all genes (by row). The per-sample weights are the distinct column values.
}
\author{Gordon Smyth}
\seealso{
This function uses output from \code{\link{glmQLFit}} and is analogous to \code{\link[limma]{arrayWeights}} in the limma package.
The \code{\link[limma]{voomWithQualityWeights}} and \code{\link[limma]{voomLmFit}} functions also estimate sample weights for count data.
}
\examples{
# Sample 2 has higher dispersion and lower weight
y <- matrix(rpois(1000*5,lambda=8),1000,5)
y[,2] <- 2 * y[,2]
fit <- glmQLFit(y)
w <- sampleWeights(fit)
barplot(as.matrix(w)[1,], main="Sample weights")
}
