\name{glmQLFTest}
\alias{glmQLFTest}

\title{Quasi-Likelihood F-Tests}

\description{Conduct genewise quasi F-tests for a given coefficient or coefficient contrast.}

\usage{
glmQLFTest(glmfit, coef = ncol(glmfit$design), contrast = NULL, fc=1, lfc=NULL,
        poisson.bound = TRUE, upshot=TRUE)
}

\arguments{
\item{glmfit}{a \code{DGEGLM} object, usually output from \code{glmQLFit}.}

\item{coef}{integer or character index vector indicating which coefficients of the linear model are to be tested equal to zero. Ignored if \code{contrast} is not \code{NULL}.}

\item{contrast}{numeric vector or matrix specifying one or more contrasts of the linear model coefficients to be tested equal to zero.}

\item{fc}{scalar or numeric vector of minimum fold-change below which changes are not considered scientifically meaningful.}

\item{lfc}{scalar or numeric vector of minimum log2-fold-change below which changes are not considered scientifically meaningful. Defaults to \code{log2(fc)}. If specified then takes precedence over \code{fc}.}

\item{poisson.bound}{logical, if \code{TRUE} then the p-value returned will never be less than would be obtained for a likelihood ratio test with NB-dispersion equal to 0 and quasi-dispersion equal to 1.
This argument is no longer relevant for the edgeR v4 QL pipeline and will be ignored if \code{glmfit} was created by \code{glmQLFit} with \code{legacy=FALSE}.}

\item{upshot}{logical. If \code{TRUE}, then a more powerful test is conducted assuming that the true logFC are uniformly distributed over the null interval from \code{-lfc} to \code{lfc}. 
If \code{FALSE}, then the standard TREAT p-values are computed assuming the worst-case scenario of true logFC equal to \code{lfc}.}
}

\details{
\code{glmQLFTest} is typically used after \code{\link{glmQLFit}}.
The two functions implement the quasi-likelihood (QL) methods of Lund et al (2012) with some enhancements and with slightly different GLM, trend and FDR methods.
See Lun et al (2016) or Chen et al (2016) for tutorials describing the use of \code{glmQLFit} and \code{glmQLFit} as part of a complete analysis pipeline.
Another case study using \code{glmQLFit} and \code{glmQLFTest} is given in Section 4.7 of the edgeR User's Guide.

\code{glmQLFTest} is similar to \code{glmLRT} except that it replaces likelihood ratio tests with empirical Bayes quasi-likelihood F-tests.
The p-values from \code{glmQLFTest} are always greater than or equal to those that would be obtained from \code{glmLRT} using the same negative binomial dispersions.
}

\value{
An object of class \code{DGELRT} with the same components as produced by \code{\link{glmLRT}}, except that the \code{table$LR} column becomes \code{table$F} and contains quasi-likelihood F-statistics.
It also stores \code{df.total}, a numeric vector containing the denominator degrees of freedom for the F-test, equal to \code{df.prior + df.residual.zeros}.
}

\references{
  Chen Y, Chen L, Lun ATL, Baldoni PL, Smyth GK (2025).
  edgeR v4: powerful differential analysis of sequencing data with expanded functionality and improved support for small counts and larger datasets.
  \emph{Nucleic Acids Research} 53(2), gkaf018.
  \doi{10.1093/nar/gkaf018}.

Chen Y, Lun ATL, Smyth GK (2016).
From reads to genes to pathways: differential expression analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline.
\emph{F1000Research} 5, 1438.
\doi{10.12688/f1000research.8987.2}.

Lun ATL, Chen Y, Smyth GK (2016).
It's DE-licious: a recipe for differential expression analyses of RNA-seq experiments using quasi-likelihood methods in edgeR.
\emph{Methods in Molecular Biology} 1418, 391-416.
\doi{10.1007/978-1-4939-3578-9_19}.
See also the preprint \url{https://gksmyth.github.io/pubs/QLedgeRPreprint.pdf}.

Lund SP, Nettleton D, McCarthy DJ, Smyth GK (2012).
Detecting differential expression in RNA-sequence data using quasi-likelihood with shrunken dispersion estimates.
\emph{Statistical Applications in Genetics and Molecular Biology} 11(5), Article 8.
\doi{10.1515/1544-6115.1826}.
See also the preprint \url{https://gksmyth.github.io/pubs/QuasiSeqPreprint.pdf}.

Lun ATL, Smyth GK (2017).
No counts, no variance: allowing for loss of degrees of freedom when assessing biological variability from RNA-seq data.
\emph{Statistical Applications in Genetics and Molecular Biology} 16(2), 83-93.
\doi{10.1515/sagmb-2017-0010}.

Phipson B, Lee S, Majewski IJ, Alexander WS, Smyth GK (2016).
Robust hyperparameter estimation protects against hypervariable genes and improves power to detect differential expression.
\emph{Annals of Applied Statistics} 10, 946-963.
\doi{10.1214/16-AOAS920}.
}

\author{Yunshun Chen, Aaron Lun, Davis McCarthy, Lizhong Chen and Gordon Smyth}

\seealso{
\code{\link{glmQLFit}} provides estimated GLMs to \code{glmQLFTest}.

\code{\link{topTags}} displays results from \code{glmQLFTest}.

\code{\link{plotQLDisp}} can be used to visualize the distribution of quasi-dispersions after EB shrinkage from \code{glmQLFit}.

The \code{QuasiSeq} package gives an alternative implementation of the Lund et al (2012) methods.
}

\examples{
# See glmQLFit for an example using glmQLFTest
}

\concept{Model fit}
\concept{Differential expression}

