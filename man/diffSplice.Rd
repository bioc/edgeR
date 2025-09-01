\name{diffSplice.DGEGLM}
\alias{diffSplice.DGEGLM}
\title{Test for Differential Transcript Usage (with edgeR Fitted Model Object)}
\description{Given an glm fit at the transcript (or exon) level, test for differential transcript (or exon) usage within genes between experimental conditions.
More generally, test for differential usage within genes of any set of splice-events or isoform-identifying features.}
\usage{
\method{diffSplice}{DGEGLM}(fit, coef=ncol(fit$design), contrast=NULL, geneid, exonid = NULL,
     robust = NULL, nexons.approx = 10, verbose = TRUE, \dots)
}
\arguments{
  \item{fit}{an \code{DGEGLM} fitted model object produced by \code{glmQLFit}. Rows should correspond to transcripts, for a DTU analysis, or to exons and exon-exon junctions for a DEU analysis.}
  \item{coef}{subscript indicating which coefficient of the generalized linear model is to be tested for differential usage. Defaults to the last coefficient. Can be an integer or a coefficient name.}
  \item{contrast}{numeric vector specifying the contrast of the linear model coefficients to be tested for differential usage. Length must equal to the number of columns of \code{design}. If specified, then takes precedence over \code{coef}.}
  \item{geneid}{gene identifiers. Either a vector of length \code{nrow(fit)} or the name of the column of \code{fit$genes} containing the gene identifiers. Rows with the same ID are assumed to belong to the same gene.}
  \item{exonid}{exon identifiers. Either a vector of length \code{nrow(fit)} or the name of the column of \code{fit$genes} containing the exon identifiers.}
  \item{robust}{logical, should the estimation of the empirical Bayes prior parameters be robustified against outlier sample variances? By default, the same setting will be used as for the \code{glmQLFit} call used to create \code{fit}.}
  \item{nexons.approx}{exact test statistics are computed for all genes with up to this number of transcripts (or exons). For genes with more exons, a more computationally fast approximation is used.}
  \item{verbose}{logical, if \code{TRUE} some diagnostic information about the number of genes and exons is output.}
  \item{\dots}{other arguments are not currently used.}
}

\value{
An object of class \code{MArrayLM} containing both exon level and gene level tests.
Results are sorted by geneid and by exonid within gene.
  \item{coefficients}{a single-column numeric matrix containing differential usage coefficients for the coefficient or contrast specified by \code{coef} or \code{contrast}. The coefficients give the difference between the log-fold-change for that genomic feature versus the average log-fold-change for all other features for the same gene.}
  \item{t}{a single-column numeric matrix of quasi t-statistics.}
  \item{p.value}{a single-column numeric vector of p-values corresponding to the t-statistics.}
  \item{genes}{data.frame of annotation for the low-level features.}
  \item{genecolname}{character string giving the name of the column of \code{genes} containing gene IDs.}
  \item{gene.F}{single-column numeric matrix of quasi F-statistics, one row for each gene.}
  \item{gene.F.p.value}{single-column numeric matrix of p-values corresponding to \code{gene.F}}
  \item{gene.simes.p.value}{single-column numeric matrix of Simes adjusted p-values, one row for each gene.}
  \item{gene.bonferroni.p.value}{single-column numeric matrix of Bonferroni adjusted p-values, one row for each gene.}
  \item{gene.genes}{data.frame of gene annotation.}
}

\details{
\code{diffSplice} is an S3 generic function defined in the limma package.
This help page describes the \code{diffSplice} method for \code{DGEGLM} fitted model objects as produced by \code{glmQLFit} in edgeR.
The output object can be further explored using the \code{topSplice} and \code{plotSplice} functions of the limma package.

The function tests for differential usage of the row-wise isoform features contained in \code{fit} for each gene across the comparison specified by the model coefficient or contrast.
The isoform features can be transcripts for a differential transcript usage (DTU) analysis, or can be a combination of exons and exon-exon junctions for a differential exon usage (DEU) analysis.

Testing for differential usage is equivalent to testing whether the log-fold-changes in the \code{fit} differ between exons for the same gene.
Two different tests are provided, one at the gene level and one at the transcript (exon) level.
In both cases, the tests are conducted using quasi-F test statistics based on deviance differences.
The gene-level test is a quasi-F-test for differences between the log-fold-changes, equivalent to testing for interaction between the contrast comparison and the transcripts for each gene.
This produces one F-statistic for each gene, where the numerator df is one less than the number of exons and the denominator is the residual df pooled over all exons for that gene.
The null hypothesis is that all the exons have the same log-fold-change and the alternative is that the log-fold-changes are not all equal.
The exon-level tests compare each exon to the other exons for the same gene.
Each exon-level test compares the model in which all the exons have the same log-fold-change vs the alternative model that just the specified exon has a different log-fold-change.
The null hypothesis is that the log-fold-change for the specified exon is equal to the consensus log-fold-change for the other exons, and the alternative is that it is not.
This produces one quasi F-statistic for each exon, where the numerator df is equal to one and the denominator is the same as for the gene-level test.
The exon-level F-statistics are then converted into t-tests by taking square-roots and multiplying by the sign of log-fold-change difference.

For genes with more than \code{nexons.approx} exons, a fast approximation is used for the exon-level test, assuming that the consensus log-fold-change for the other exons should be close to the consensus log-fold-change for all exons for that gene.
The null hypothesis for these genes is that the log-fold-change for the specified exon is equal to the consensus log-fold-change for all other exons, and the alternative is that it is not.
The approximation should be very slightly conservative.

The \code{fit} object should be created using \code{glmQLFit}.
\code{diffSplice} will automatically detect whether \code{glmQLFit} was run with edgeR v4 bias-adjusted deviances (\code{legacy=FALSE}, recommended) or with ordinary deviances (\code{legacy=TRUE}).

The exon-level tests are converted into genewise tests by adjusting the p-values for the same gene by Simes method.
Alternatively, the exon-level tests are also converted into genewise tests by adjusting the smallest p-value for each gene by Bonferroni's method.

This function can be used on transcript level RNA-seq counts from Salmon or kallisto, after using catchSalmon() or catchKallisto() and glmQLFit(), as described by Baldoni et al (2025).
It can also be used on equivalence-class counts from Salmon or kallisto, after using glmQLFit(), as described by Cmero et al (2019).
It can also be used on exon-level read counts.
}

\note{
\code{\link{diffSpliceDGE}} is a legacy function for the same purpose as \code{diffSplice.DGEGLM}, and \code{diffSplice.DGEGLM} is intended to replace it.
\code{diffSplice.DGEGLM} has better statistical power and receiver operating curve (ROC), while still controlling the FDR.
In terms of speed, the new function is slightly slower but the difference should not be enough to have a material effect on analyses.
}

\seealso{
\code{\link{diffSpliceDGE}} is an older function with a similar purpose.
See \code{\link{diffSplice}}, \code{\link{topSplice}}, and \code{\link{plotSplice}} in the limma package.
}

\author{Lizhong Chen, Yunshun Chen and Gordon Smyth}

\references{
Baldoni PL, Chen L, Li M, Chen Y, Smyth GK (2025).
Dividing out quantification uncertainty enables assessment of differential transcript usage with limma and edgeR.
\emph{bioRxiv}
\doi{10.1101/2025.04.07.647659}.

Cmero M, Davidson NM, Oshlack A (2019).
Using equivalence class counts for fast and accurate testing of differential transcript usage.
\emph{F1000Research} 8, 265.
\doi{10.12688/f1000research.18276.2}.
}

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
\concept{Differential transcript usage}
\concept{Differential usage}
