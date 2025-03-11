\name{diffSplice.DGEGLM}
\alias{diffSplice.DGEGLM}
\title{Test for Differential Splicing Using edgeR Fitted Model}
\description{Given a linear model fit at the exon (transcript) level, test for differences in exon (transcript) usage between experimental conditions.}
\usage{
\method{diffSplice}{DGEGLM}(fit, coef=ncol(fit$design), contrast=NULL, geneid, exonid = NULL,
     robust = FALSE, nexons.approx = 10, verbose = TRUE, \dots)
}
\arguments{
  \item{fit}{an \code{DGEGLM} fitted model object produced by \code{glmQLFit}. Rows should correspond to low-level genomic features such as exons, exon-exon junctions or transcripts.}
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
This function tests for differential usage of the genomic features for each gene across the comparison specified by the model coefficient or contrast.
The genomic features are usually transcripts, exons, or exon-exon junctions.

Testing for differential exon (transcript) usage is equivalent to testing whether the log-fold-changes in the \code{fit} differ between exons for the same gene.
Two different tests are provided, one at the gene level and one at the exon (transcript) level.
In both cases, the tests are conducted using quasi-F test statistics based on deviance differences.
The gene-level test is a quasi-F-test for differences between the log-fold-changes, equivalent to testing for interaction between the contrast comparison and the transcripts for each gene.
This produces one F-statistic for each gene, where the numerator df is one less than the number of exons and the denominator is the residual df pooled over all exons for that gene.
The null hypothesis is that all the exons have the same log-fold-change and the alternative is that the log-fold-changes are not all equal.
The exon-level tests compare each exon to the other exons for the same gene.
Each exon-level test compares the model in which all the exons have the same log-fold-change vs the alternative model that just the specified exon has a different log-fold-change.
The null hypothesis is that the log-fold-change for the specified exon is equal to the consensus log-fold-change for the other exons, and the alternative is that it is not.
This produces one F-statistic for each exon, where the numerator df is equal to one and the denominator is the same as for the gene-level test.
The exon-level F-statistics are then converted into t-tests by taking square-roots and multiplying by the sign of log-fold-change difference.

For genes with more than \code{nexons.approx} exons, a fast approximation is used for the exon-level test by assuming that the consensus log-fold-change for the other exons is the same as the consensus log-fold-change for all exons for that gene.
The approximation should be very slightly conservative.

The \code{fit} object should be created using \code{glmQLFit}.
\code{diffSplice} will automatically detect whether \code{glmQLFit} was run with edgeR v4 bias-adjusted deviances (\code{legacy=FALSE}, recommended) or with ordinary deviances (\code{legacy=TRUE}).

The exon-level tests are converted into genewise tests by adjusting the p-values for the same gene by Simes method.
Alternatively, the exon-level tests are also converted into genewise tests by adjusting the smallest p-value for each gene by Bonferroni's method.

The \code{diffSplice} function is an S3 generic and will accept either a linear model fit from limma or a glm QL fit from edgeR.
The output object is the same in either case.
}

\seealso{
\code{\link{diffSpliceDGE}} is an older function with a similar purpose.
See \code{\link{diffSplice}}, \code{\link{topSplice}}, and \code{\link{plotSplice}} in the limma package.
}

\author{Lizhong Chen, Yunshun Chen and Gordon Smyth}

\examples{
# Gene and exon annotation
Gene <- paste("Gene", 1:100, sep="")
Gene <- rep(Gene, each=10)
Exon <- paste("Ex", 1:10, sep="")
Gene.Exon <- paste(Gene, Exon, sep=".")
genes <- data.frame(GeneID=Gene, Gene.Exon=Gene.Exon)

# Two groups with n=3 replicates in each group
group <- factor(rep(1:2, each=3))
design <- model.matrix(~group)

# Generate exon counts.
# Knock-out the first exon of Gene1 by 90%
mu <- matrix(100, nrow=1000, ncol=6)
mu[1,4:6] <- 10
counts <- matrix(rnbinom(6000,mu=mu,size=20),1000,6)
y <- DGEList(counts=counts, lib.size=rep(1e6,6), genes=genes)

# Exon level fit
fit <- glmQLFit(y, design)

# Differential usage analysis
ds <- diffSplice(fit, geneid="GeneID")
topSplice(ds)
plotSplice(ds)
}

\keyword{rna-seq}
\concept{Differential exon usage}
