\name{normLibSizes}
\alias{normLibSizes}
\alias{normLibSizes.DGEList}
\alias{normLibSizes.SummarizedExperiment}
\alias{normLibSizes.default}
\alias{calcNormFactors}
\alias{calcNormFactors.DGEList}
\alias{calcNormFactors.SummarizedExperiment}
\alias{calcNormFactors.default}

\title{Normalize Library Sizes}
\description{
Calculate scaling factors to convert the raw library sizes for a set of sequenced samples into normalized effective library sizes.
}
\usage{
\method{normLibSizes}{DGEList}(object,
                method = c("TMM", "TMMwsp", "RLE", "upperquartile", "none"),
                refColumn = NULL, logratioTrim = .3, sumTrim = 0.05, doWeighting = TRUE,
                Acutoff = -1e10, p = 0.75, \dots)
\method{normLibSizes}{SummarizedExperiment}(object,
                method = c("TMM","TMMwsp","RLE","upperquartile","none"),
                refColumn = NULL, logratioTrim = .3, sumTrim = 0.05, doWeighting = TRUE,
                Acutoff = -1e10, p = 0.75, \dots)
\method{normLibSizes}{default}(object, lib.size = NULL,
                method = c("TMM","TMMwsp","RLE","upperquartile","none"),
                refColumn = NULL, logratioTrim = .3, sumTrim = 0.05, doWeighting = TRUE,
                Acutoff = -1e10, p = 0.75, \dots)

\method{calcNormFactors}{DGEList}(object,
                method = c("TMM","TMMwsp","RLE","upperquartile","none"),
                refColumn = NULL, logratioTrim = .3, sumTrim = 0.05, doWeighting = TRUE,
                Acutoff = -1e10, p = 0.75, \dots)
\method{calcNormFactors}{SummarizedExperiment}(object,
                method = c("TMM","TMMwsp","RLE","upperquartile","none"),
                refColumn = NULL, logratioTrim = .3, sumTrim = 0.05, doWeighting = TRUE,
                Acutoff = -1e10, p = 0.75, \dots)
\method{calcNormFactors}{default}(object, lib.size = NULL,
                method = c("TMM","TMMwsp","RLE","upperquartile","none"),
                refColumn = NULL, logratioTrim = .3, sumTrim = 0.05, doWeighting = TRUE,
                Acutoff = -1e10, p = 0.75, \dots)
}

\arguments{
  \item{object}{a matrix of raw read counts or a \code{DGEList} object or a \code{SummarizedExperiment} object.}
  \item{lib.size}{numeric vector of library sizes corresponding to the columns of the matrix \code{object}.}
  \item{method}{normalization method to be used.}
  \item{refColumn}{column to use as reference for \code{method="TMM"}. Can be a column number or a numeric vector of length \code{nrow(object)}.}
  \item{logratioTrim}{the fraction (0 to 0.5) of observations to be trimmed from each tail of the distribution of log-ratios (M-values) before computing the mean. Used by \code{method="TMM"} for each pair of samples.}
  \item{sumTrim}{the fraction (0 to 0.5) of observations to be trimmed from each tail of the distribution of A-values before computing the mean. Used by \code{method="TMM"} for each pair of samples.}
  \item{doWeighting}{logical, whether to use (asymptotic binomial precision) weights when computing the mean M-values. Used by \code{method="TMM"} for each pair of samples.}
  \item{Acutoff}{minimum cutoff applied to A-values. Count pairs with lower A-values are ignored. Used by \code{method="TMM"} for each pair of samples.}
  \item{p}{numeric value between 0 and 1 specifying which quantile of the counts should be used by \code{method="upperquartile"}.}
  \item{\dots}{other arguments are not currently used.}
}

\details{
This function computes scaling factors to convert observed library sizes into normalized library sizes, also called "effective library sizes".
The effective library sizes for use in downstream analysis are \code{lib.size * norm.factors} where \code{lib.size} contains the original library sizes and \code{norm.factors} is the vector of scaling factors computed by this function.

The TMM method implements the trimmed mean of M-values method proposed by Robinson and Oshlack (2010).
By default, the M-values are weighted according to inverse variances, as computed by the delta method for logarithms of binomial random variables.
If \code{refColumn} is unspecified, then the column whose count-per-million upper quartile is closest to the mean upper quartile is set as the reference library.

The TMMwsp method stands for "TMM with singleton pairing".
This is a variant of TMM that is intended to have more stable performance when the counts have a high proportion of zeros.
In the TMM method, genes that have zero count in either library are ignored when comparing pairs of libraries.
In the TMMwsp method, the positive counts from such genes are reused to increase the number of features by which the libraries are compared.
The singleton positive counts are paired up between the libraries in decreasing order of size and then a slightly modified TMM method is applied to the re-ordered libraries.
If \code{refColumn} is unspecified, then the column with largest sum of square-root counts is used as the reference library.

RLE is the scaling factor method proposed by Anders and Huber (2010).
We call it "relative log expression", as median library is calculated from the geometric mean of all columns and the median ratio of each sample to the median library is taken as the scale factor.

The upperquartile method is the upper-quartile normalization method of Bullard et al (2010), in which the scale factors are calculated from the 75\% quantile of the counts for each library, after removing genes that are zero in all libraries.
The idea is generalized here to allow normalization by any quantile of the count distributions.

If \code{method="none"}, then the normalization factors are set to 1.

For symmetry, normalization factors are adjusted to multiply to 1.
Rows of \code{object} that have zero counts for all columns are removed before normalization factors are computed.
The number of such rows does not affect the estimated normalization factors.

If \code{object} is a \code{SummarizedExperiment} object, then it is converted to a \code{DGEList} using \code{SE2DGEList} and the \code{DGEList} method applied.
}

\value{
If \code{object} is a \code{matrix}, then the output is a vector with length \code{ncol(object)} giving the library scaling factors.

If \code{object} is a \code{DGEList} or \code{SummarizedExperiment} object, then the output is a \code{DGEList} the same as input with the library scaling factors stored as \code{object$samples$norm.factors}.
}

\note{
\code{normLibSizes} is the new name for \code{calcNormFactors}.
The two functions are equivalent but \code{calcNormFactors} will eventually be retired.
}

\author{Mark Robinson, Gordon Smyth, Yunshun Chen.}

\references{
Anders, S, Huber, W (2010).
Differential expression analysis for sequence count data
\emph{Genome Biology} 11, R106.

Bullard JH, Purdom E, Hansen KD, Dudoit S. (2010)
Evaluation of statistical methods for normalization and differential expression in mRNA-Seq experiments.
\emph{BMC Bioinformatics} 11, 94.

Robinson MD, Oshlack A (2010).
A scaling normalization method for differential expression analysis of RNA-seq data.
\emph{Genome Biology} 11, R25.
}

\seealso{
\code{\link{getNormLibSizes}}, \code{\link{SE2DGEList}}.
}

\examples{
y <- matrix( rpois(1000, lambda=5), nrow=200 )
normLibSizes(y)

# The TMM and TMMwsp methods give zero weight to the genes with the largest fold-changes:
y <- cbind(1,c(1,1,1,1,1,1,1,1,1,100))
normLibSizes(y, lib.size=c(1e6,1e6))

# normLibSizes makes the fold-changes for the majority of genes as close to 1 as possible:
# In this example, computing CPMs with raw library sizes makes the two samples look quite different.
dge <- DGEList(counts=y)
cpm(dge)
# By contrast, normalizing the library sizes makes most of the CPMs equal in the two samples:
dge <- normLibSizes(dge)
cpm(dge)
getNormLibSizes(dge)
}

\concept{Normalization}
