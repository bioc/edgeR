\name{diffSplice.DGEGLM}
\alias{diffSplice.DGEGLM}
\title{Test for Differential Splicing Using edgeR Fitted Model}
\description{Given a linear model fit at the exon (transcript) level, test for differences in exon (transcript) usage between experimental conditions.}
\usage{
\method{diffSplice}{DGEGLM}(fit, coef=ncol(fit$design), contrast=NULL, geneid, exonid = NULL, robust = FALSE, nexons.approx = 10, verbose = TRUE, \dots)
}
\arguments{
  \item{fit}{an \code{DGEGLM} fitted model object produced by \code{glmQLFit}. Rows should correspond to exons.}
  \item{coef}{integer indicating which coefficient of the generalized linear model is to be tested for differential exon usage. Defaults to the last coefficient.}
  \item{contrast}{numeric vector specifying the contrast of the linear model coefficients to be tested for differential exon usage. Length must equal to the number of columns of \code{design}. If specified, then takes precedence over \code{coef}.}
  \item{geneid}{gene identifiers. Either a vector of length \code{nrow(fit)} or the name of the column of \code{fit$genes} containing the gene identifiers. Rows with the same ID are assumed to belong to the same gene.}
  \item{exonid}{exon identifiers. Either a vector of length \code{nrow(fit)} or the name of the column of \code{fit$genes} containing the exon identifiers.}
  \item{robust}{logical, should the estimation of the empirical Bayes prior parameters be robustified against outlier sample variances?}
  \item{nexons.approx}{exact test statistics are computed for all genes with up to this number of exons. For genes with more exons, a more computationally fast approximation is used.}
  \item{verbose}{logical, if \code{TRUE} some diagnostic information about the number of genes and exons is output.}
  \item{\dots}{other arguments are not currently used.}
}
\value{
An object of class \code{MArrayLM} containing both exon level and gene level tests.
Results are sorted by geneid and by exonid within gene.
  \item{coefficients}{numeric matrix of coefficients of same dimensions as \code{fit}. Each coefficient is the difference between the log-fold-change for that exon versus the average log-fold-change for all other exons for the same gene.}
  \item{t}{numeric matrix of moderated t-statistics, of same dimensions as \code{fit}.}
  \item{p.value}{numeric vector of p-values corresponding to the t-statistics}
  \item{genes}{data.frame of exon annotation}
  \item{genecolname}{character string giving the name of the column of \code{genes} containing gene IDs}
  \item{gene.F}{numeric matrix of moderated F-statistics, one row for each gene.}
  \item{gene.F.p.value}{numeric matrix of p-values corresponding to \code{gene.F}}
  \item{gene.simes.p.value}{numeric matrix of Simes adjusted p-values, one row for each gene.}
  \item{gene.bonferroni.p.value}{numeric matrix of Bonferroni adjusted p-values, one row for each gene.}
  \item{gene.genes}{data.frame of gene annotation.}
}
\details{
This function tests for differential usage of the genomic features for each gene and for each column of \code{fit}.
The genomic features are usually transcripts, or exons, or exon-exon junctions.

Testing for differential exon usage is equivalent to testing whether the log-fold-changes in the \code{fit} differ between exons for the same gene.
Two different tests are provided.
The first is an F-test for differences between the log-fold-changes.
The other is a series of t-tests in which each exon is compared to the average of all other exons for the same gene.
The exon-level t-tests are converted into genewise tests by adjusting the p-values for the same gene by Simes method.
Alternatively, the exon-level t-tests are also converted into genewise tests by adjusting the smallest p-value for each gene by Bonferroni's method.
The minimum adjusted p-value is then used for each gene.

The \code{diffSplice} function is S3 generic and will accept a linear model fit from limma or a glm QL fit from edgeR.
The output object is the same in either case.
}
\seealso{
\code{\link{diffSplice}}, \code{\link{topSplice}}, \code{\link{plotSplice}}
}
\author{Gordon Smyth and Charity Law}

\examples{
\dontrun{
v <- voom(dge,design)
fit <- lmFit(v,design)
ex <- diffSplice(fit,geneid="EntrezID")
topSplice(ex)
plotSplice(ex)
}
}

\keyword{rna-seq}
