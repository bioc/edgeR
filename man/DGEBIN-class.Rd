\name{DGEBIN-class}
\docType{class}
\alias{DGEBIN-class}
\alias{show,DGEBIN-method}
\title{Digital Gene Expression Binomial Model Results - Class}

\description{
A list-based S4 class for storing results of a binomial generalized linear model fit to each gene in a DGE dataset, as produced by \code{\link{binFit}} or \code{\link{binQLFit}}.
It is a sibling of \code{\link{DGEGLM-class}} (the negative-binomial GLM fit): a binomial fit is deliberately \emph{not} a \code{DGEGLM}, so that binomial-specific tests such as \code{\link{binQLFTest}} accept only a \code{DGEBIN} and the negative-binomial tests (\code{\link{glmLRT}}, \code{\link{glmQLFTest}}, \code{\link{glmTreat}}) do not.
}

\section{List Components}{
For objects of this class, rows correspond to genomic features and columns to coefficients in the linear model.

Objects of this class contain the following list components:
\describe{
  \item{\code{coefficients}:}{matrix of estimated coefficients from the binomial GLM fits, on the natural log scale.}
  \item{\code{unshrunk.coefficients}:}{matrix of coefficients from the unshrunk fits; present only when \code{prior.count>0}.}
  \item{\code{fitted.values}:}{matrix of fitted (expected) proportions from the binomial GLM for each gene.}
  \item{\code{deviance}:}{vector giving the deviance from the model fit to each gene.}
  \item{\code{df.residual}:}{vector of residual degrees of freedom for the model fit to each gene.}
  \item{\code{counts}:}{the first matrix of counts (successes) as input.}
  \item{\code{counts2}:}{the second matrix of counts (failures) as input.}
  \item{\code{design}:}{design matrix for the full model.}
  \item{\code{weights}:}{matrix of prior weights used in the GLM fitting.}
  \item{\code{nthreads}:}{integer giving the number of OpenMP threads used for the fit.}
  \item{\code{AveLogCPM}:}{numeric vector giving average log2 counts per million of the total coverage for each gene.}
  \item{\code{samples}, \code{genes}:}{sample and gene annotation data frames, attached when the input is a \code{\link{PCList-class}} object.}
}
Objects produced by \code{\link{binQLFit}} additionally carry the leverage-adjusted and empirical Bayes quasi-dispersion components: \code{df.residual.adj}, \code{deviance.adj}, \code{df.prior}, \code{s2.post}, \code{s2.prior} and \code{average.ql.dispersion}.
}

\section{Methods}{
This class inherits directly from class \code{list} so any operation appropriate for lists will work on objects of this class.

The dimensions, row names and column names of a \code{DGEBIN} object are defined by those of the fit, see \code{\link{dim.DGEBIN}} or \code{\link{dimnames.DGEBIN}}.

\code{DGEBIN} objects can be subsetted by row, see \code{\link{subsetting}}.

\code{DGEBIN} objects also have a \code{show} method so that printing produces a compact summary of their contents.
}

\author{Lizhong Chen, Gordon Smyth}

\seealso{
\code{\link{binFit}} and \code{\link{binQLFit}} produce \code{DGEBIN} objects; \code{\link{binQLFTest}} tests them.

Other classes defined in edgeR are \code{\link{DGEGLM-class}}, \code{\link{DGEList-class}}, \code{\link{DGEExact-class}}, \code{\link{DGELRT-class}}, \code{\link{TopTags-class}}.
}

\concept{edgeR classes}
