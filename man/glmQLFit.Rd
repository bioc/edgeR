\name{glmQLFit}
\alias{glmQLFit}
\alias{glmQLFit.DGEList}
\alias{glmQLFit.SummarizedExperiment}
\alias{glmQLFit.default}

\title{Genewise Negative Binomial Generalized Linear Models with Quasi-Dispersion Estimation}

\description{Fit a negative binomial generalized log-linear model to the read counts for each gene.
Estimate the genewise quasi-dispersions with empirical Bayes moderation.}

\usage{
\method{glmQLFit}{default}(y, design = NULL, dispersion = NULL, offset = NULL, lib.size = NULL,
        weights = NULL, abundance.trend = TRUE, AveLogCPM = NULL,
        covariate.trend = NULL, robust = FALSE, winsor.tail.p = c(0.05, 0.1),
        legacy = FALSE, nthreads = 1L, \dots)
\method{glmQLFit}{DGEList}(y, design = NULL, dispersion = NULL, abundance.trend = TRUE,
        robust = FALSE, winsor.tail.p = c(0.05, 0.1),
        legacy = FALSE, \dots)
\method{glmQLFit}{SummarizedExperiment}(y, design = NULL, dispersion = NULL, abundance.trend = TRUE,
        robust = FALSE, winsor.tail.p = c(0.05, 0.1),
        legacy = FALSE, \dots)
}

\arguments{
\item{y}{a matrix of counts or a \code{DGEList} object or a \code{SummarizedExperiment} containing counts. Rows represent genes and columns represent samples.}

\item{design}{
  numeric matrix giving the design matrix for the genewise linear models.
  Must be of full column rank.
  Defaults to a single column of ones, equivalent to treating the columns as replicate libraries.
}

\item{dispersion}{
  negative binomial dispersions. Can be a single value, or a vector of length \code{nrow(y)}, or a matrix of the same size as \code{y}.
  In the new edgeR v4 pipeline (\code{legacy=FALSE}), is usually left \code{NULL}. Only preset \code{dispersion} if you know what you are doing! In this pipeline, any dispersion values greater than 4 will be capped at 4.
  In the older edgeR v3 pipeline (\code{legacy=TRUE}), will be extracted from the \code{DGEList} object \code{y}. The trended dispersions will be extracted if present, otherwise common dispersion, or will be set to 0.05 if neither is present.
} 

\item{offset}{offsets for the log-linear models containing log effective library sizes.  Can be a single value, or a vector of length \code{ncol(y)}, or a matrix of the same size as \code{y}. If \code{NULL}, then will be computed by \code{getOffset(y)}.}

\item{lib.size}{numeric vector of length \code{ncol(y)} giving library sizes. Only used if \code{offset=NULL}, in which case \code{offset} is set to \code{log(lib.size)}. Defaults to \code{colSums(y)}.}

\item{weights}{positive prior weights for the GLM fits. Can be a single value, or a vector of length \code{ncol(y)}, or a matrix of the same size as \code{y}. If \code{NULL}, will be set to unity for all observations.}

\item{abundance.trend}{logical, whether to allow an abundance trend for the quasi-dispersion prior.}

\item{AveLogCPM}{numeric vector giving average log2-counts per million for each row of \code{y}. If \code{NULL}, then will be computed by \code{aveLogCPM(y)}.}

\item{covariate.trend}{numeric vector of length \code{nrow(y)}. If non-NULL, then will be used instead of \code{AveLogCPM} as the covariate for the trended quasi-dispersion prior.}

\item{robust}{logical, whether to estimate the prior distribution for the quasi-dispersions robustly.}

\item{winsor.tail.p}{numeric vector of length 2 giving proportion to trim (Winsorize) from lower and upper tail of the distribution of genewise deviances when estimating the hyperparameters of the quasi-dispersion prior. Used as input to \code{squeezeVar} when \code{robust=TRUE}.} 

\item{legacy}{logical, if \code{TRUE} then produce legacy results as for Bioconductor 3.16 and earlier. If \code{FALSE}, then use the new quasi-dispersion estimation method with adjusted deviances.}

\item{nthreads}{integer giving the number of OpenMP threads to use for the genewise computations. Defaults to 1 (serial).}

\item{\dots}{other arguments are passed to \code{\link{glmFit}}.}
}

\details{
\code{glmQLFit} and \code{glmQLFTest} implement the quasi-likelihood (QL) methods of Lund et al (2012) with some enhancements and with slightly different GLM, trend and FDR methods.
See Lun et al (2016) or Chen et al (2016) for tutorials describing the use of \code{glmQLFit} and \code{glmQLFit} as part of a complete analysis pipeline.
Another case study using \code{glmQLFit} and \code{glmQLFTest} is given in Section 4.7 of the edgeR User's Guide.

\code{glmQLFit} is similar to \code{glmFit} except that it also estimates a quasi-dispersion (quasi-dispersion) for each gene.
It calls the limma function \code{\link[limma]{squeezeVar}} to conduct empirical Bayes moderation of the genewise quasi-dispersions.
If \code{robust=TRUE}, then the robust hyperparameter estimation features of \code{squeezeVar} are used (Phipson et al, 2016).
If \code{abundance.trend=TRUE}, then a prior trend is estimated based on the average logCPMs.
If \code{covariate.trend} is not \code{NULL}, then a prior trend is estimated using the \code{covariate.trend} values as predictors.

\code{glmQLFit} gives special attention to handling of small counts and zero counts.
When \code{legacy=TRUE}, the function uses the method of Lun and Smyth (2017) to adjust the residual degrees of freedom when fitted values of zero provide no useful residual degrees of freedom for estimating the quasi-dispersion.
The usual residual degrees of freedom are returned as \code{df.residual} while the adjusted residual degrees of freedom are returned as \code{df.residuals.zeros}.

If \code{legacy=FALSE}, then a more comprehensive adjustment for small counts is used.
The new method adjusts both the residual deviances and the residual degrees of freedom to improve the accuracy of the quasi-dispersion estimates even for very small counts (Chen et al 2024).
With this new method, the residual deviance is no longer equal to the usual residual deviance from generalized linear model theory and the residual degrees of freedom are no longer integers.
With the new method, the \code{glmQLFTest} function should give good error rate control even for data with many small counts.
The \code{legacy=FALSE} method was introduced in edgeR 4.0.0 with the Bioconductor 3.18 release.
Setting \code{legacy=TRUE} will reproduce earlier results as for edgeR 3.8.0 and Bioconductor 3.16.

\code{glmQLFit} requires the NB dispersion to be pre-specified when \code{legacy=TRUE} but not when \code{legacy=FALSE}.
In the new pipeline, the NB dispersion will be automatically estimated by \code{glmQLFit} whenever \code{dispersion=NULL}.
In that case, the NB-dispersion will be estimated as the common dispersion of the most highly expressed genes by \code{AveLogCPM}.
The number of top genes used is \code{min(ceiling(p*ngenes), 10000)}, where
\code{p = chooseLowessSpan(ngenes*sqrt(df.residual),small.n=20,min.span=0.02)},
\code{ngenes} is the number of genes (rows of \code{y}) and \code{df.residual} is the residual degrees of freedom, equal to \code{ncol(y)-ncol(design)}.
}

\note{
The negative binomial dispersions \code{dispersion} supplied to \code{glmQLFit} must be based on a global model, that is, they must be either trended or common dispersions.
It is not correct to supply genewise dispersions because \code{glmQLFit} estimates gene-specific variability using the quasi-dispersions.
}

\value{
\code{glmQLFit} with \code{legacy=TRUE} produces an object of class \code{DGEGLM} with the same components as produced by \code{\link{glmFit}}, plus:
\item{df.residual.zeros}{numeric vector of effective residual degrees of freedom for each gene, after accounting for treatment groups with all zero counts.}
\item{df.prior}{numeric vector of prior degrees of freedom for the quasi-dispersions. Has length \code{nrow(y)} if \code{robust=TRUE}, otherwise length 1.}
\item{s2.prior}{numeric vector giving prior value for the quasi-dispersions. Has length \code{nrow(y)} if \code{abundance.trend=TRUE}, otherwise length 1.}
\item{s2.post}{numeric vector of posterior genewise quasi-dispersions.}

\code{glmQLFit} with \code{legacy=FALSE} produces an object of class \code{DGEGLM} with the same components as produced by \code{\link{glmFit}}, plus:
\item{df.residual.adj}{numeric vector of adjusted residual degrees of freedom for each gene.}
\item{df.prior}{numeric vector of prior degrees of freedom for the quasi-dispersions. Has length \code{nrow(y)} if \code{robust=TRUE}, otherwise length 1.}
\item{s2.prior}{numeric vector giving prior value for the quasi-dispersions. Has length \code{nrow(y)} if \code{abundance.trend=TRUE}, otherwise length 1.}
\item{s2.post}{numeric vector of posterior genewise quasi-dispersions.} 
\item{average.ql.dispersion}{average quasi-dispersion, used to scale the NB dispersion after estimating the quasi-dispersions.}
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
\url{https://gksmyth.github.io/pubs/QuasiSeqPreprint.pdf}.

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
\code{\link{glmQLFTest}} performs F-tests using the fit from \code{glmQLFit}.

\code{\link{plotQLDisp}} can be used to visualize the distribution of quasi-dispersions after EB shrinkage from \code{glmQLFit}.

The \code{QuasiSeq} package gives an alternative implementation of the Lund et al (2012) methods.
}

\examples{
nlibs <- 4
ngenes <- 1000
dispersion.true <- 1/rchisq(ngenes, df=10)
design <- model.matrix(~factor(c(1,1,2,2)))

# Generate count data
y <- rnbinom(ngenes*nlibs,mu=20,size=1/dispersion.true)
y <- matrix(y,ngenes,nlibs)
d <- DGEList(y)
d <- normLibSizes(d)

# Fit the NB GLMs with QL methods
fit <- glmQLFit(d, design)
results <- glmQLFTest(fit)
topTags(results)
fit <- glmQLFit(d, design, robust=TRUE)
results <- glmQLFTest(fit)
topTags(results)
fit <- glmQLFit(d, design, abundance.trend=FALSE)
results <- glmQLFTest(fit)
topTags(results)
}

\concept{Model fit}
\concept{Dispersion estimation}
\concept{Differential expression}

