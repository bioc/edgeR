\name{decideTests}
\alias{decideTestsDGE}
\alias{decideTests.DGEExact}
\alias{decideTests.DGELRT}
\title{Multiple Testing Across Genes and Contrasts}
\description{
Identify which genes are significantly differentially expressed from an edgeR test object containing p-values and test statistics.
}
\usage{
decideTestsDGE(object, adjust.method="BH", p.value=0.05, lfc=0)
\method{decideTests}{DGELRT}(object, adjust.method="BH", p.value=0.05, lfc=0, \dots)
}
\arguments{
  \item{object}{an object of class \code{DGEExact}, \code{DGELRT} or \code{glmQLFTest} from which p-values and log-fold-changes can be extracted.}
  \item{adjust.method}{character string specifying p-value adjustment method.
  Possible values are \code{"none"}, \code{"BH"}, \code{"fdr"}, \code{"BY"} and \code{"holm"}.
  See \code{\link[stats]{p.adjust}} for details.}
  \item{p.value}{numeric value between 0 and 1 giving the required family-wise error or false discovery rate.}
  \item{lfc}{numeric, minimum absolute log2-fold-change required.}
  \item{\dots}{other arguments are not used.}
}
\value{
An object of class \code{\link[limma:TestResults]{TestResults}}.
This is essentially a single-column integer matrix with elements \code{-1}, \code{0}
or \code{1} indicating whether each gene is classified as
significantly down-regulated, not significant or
significant up-regulated for the comparison contained in \code{object}.
To be considered significant, genes need to have adjusted p-value below \code{p.value} and log2-fold-change greater than \code{lfc}.

If \code{object} contains F-tests or LRTs for multiple contrasts, then the genes are simply classified as significant (1) or not significant.
In this case, the log2-fold-change theshold \code{lfc} has to be achieved by at least one of the contrasts for a gene to be significant.
}
\details{
This function applies a multiple testing procedure and significance level cutoff to the genewise tests contained in an edgeR test object and collates the results in a data.frame table.

The function can apply optionally apply a logFC cutoff and well as a p-value or FDR cutoff (although logFCs cutoff are not recommended, see note below).
If the statistical tests are on 1 degree of freedom, then the logFC cutoff is applied to the absolute coefficient or contrast.
If the statistical tests are on more than 1 degree of freedom, then the logFC cutoff will be satisfied if any of the coefficients or contrasts that define the test are greater than the cutoff.

\code{decideTests} and \code{decideTestDGE} are equivalent.
\code{decideTestsDGE} is currently retained for backward compatibility but will be removed in a future version of edgeR.
}
\note{
Although this function enables users to set p-value and logFC cutoffs simultaneously, this combination criterion is not recommended.
logFC cutoffs tend to favor low expressed genes and thereby reduce rather than increase biological significance.
Unless the fold changes and p-values are very highly correlated, the addition of a fold change cutoff can also increase the family-wise error rate or false discovery rate above the nominal level.
Users wanting to use fold change thresholding should considering using \code{glmTreat} instead and leaving \code{lfc} at the default value when using \code{decideTests}.
}
\seealso{
\code{\link[limma:decideTests]{decideTests}} and \code{\link[limma:TestResults]{TestResults}} in the limma package.
}
\author{Davis McCarthy, Gordon Smyth and the edgeR team}
\examples{
ngenes <- 100
x1 <- rnorm(6)
x2 <- rnorm(6)
design <- cbind(Intercept=1,x1,x2)
beta <- matrix(0,ngenes,3)
beta[,1] <- 4
beta[1:20,2] <- rnorm(20)
mu <- 2^(beta \%*\% t(design))
y <- matrix(rnbinom(ngenes*6,mu=mu,size=10),ngenes,6)
fit <- glmFit(y,design,dispersion=0.1)
lrt <- glmLRT(fit,coef=2:3)
res <- decideTests(lrt,p.value=0.1)
summary(res)
lrt <- glmLRT(fit,coef=2)
res <- decideTests(lrt,p.value=0.1)
summary(res)
}

\concept{Differential expression}
