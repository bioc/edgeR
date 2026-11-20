\name{glmFit}
\alias{glmFit}
\alias{glmFit.DGEList}
\alias{glmFit.SummarizedExperiment}
\alias{glmFit.default}

\title{Genewise Negative Binomial Generalized Linear Models}

\description{Fit a negative binomial generalized log-linear model to the read counts for each gene.}

\usage{
\method{glmFit}{default}(y, design = NULL, dispersion = NULL,
       offset = NULL, lib.size = NULL, weights = NULL,
       prior.count = 0.125, start = NULL, nthreads = 1, \dots)
\method{glmFit}{DGEList}(y, design = NULL, dispersion = NULL,
       prior.count = 0.125, start = NULL, nthreads = 1, \dots)
\method{glmFit}{SummarizedExperiment}(y, design = NULL, dispersion = NULL,
       prior.count = 0.125, start = NULL, nthreads = 1, \dots)
}

\arguments{
\item{y}{a matrix of counts or a \code{DGEList} object or a \code{SummarizedExperiment} containing counts. Rows represent genes and columns represent samples.}

\item{design}{
  numeric matrix giving the design matrix for the genewise linear models.
  Must be of full column rank.
  Defaults to a single column of ones for \code{glmFit.default} or to a oneway layout defined by \code{y$samples$group} for \code{glmFit.DGEList}.
}

\item{dispersion}{negative binomial dispersions. Can be a single value, or a vector of length \code{nrow(y)}, or a matrix of the same size as \code{y}. If \code{NULL}, will be extracted from \code{y} in this order of precedence: genewise dispersion, trended dispersions, common dispersion.}

\item{offset}{offsets for the log-linear models containing log effective library sizes.  Can be a single value, or a vector of length \code{ncol(y)}, or a matrix of the same size as \code{y}. If \code{NULL}, then will be computed by \code{getOffset(y)}.}

\item{lib.size}{numeric vector of length \code{ncol(y)} giving library sizes. Only used if \code{offset=NULL}, in which case \code{offset} is set to \code{log(lib.size)}. Defaults to \code{colSums(y)}.}

\item{weights}{positive prior weights for the GLM fits. Can be a single value, or a vector of length \code{ncol(y)}, or a matrix of the same size as \code{y}. If \code{NULL}, will be set to unity for all observations.}

\item{prior.count}{average prior count to be added to each observation to shrink the estimated log-fold-changes towards zero.}

\item{start}{optional numeric matrix of initial estimates for the linear model coefficients.}

\item{nthreads}{number of OpenMP threads to use.}

\item{\dots}{other arguments are not currently used.}
}

\value{
An object of class \code{DGEGLM} containing components \code{counts}, \code{samples}, \code{genes} and \code{abundance} from \code{y} plus the following new components:
	\item{design}{design matrix as input.}
	\item{weights}{matrix of weights as input.}
	\item{df.residual}{numeric vector of residual degrees of freedom, one for each gene.}
	\item{offset}{numeric matrix of linear model offsets.}
	\item{dispersion}{vector of dispersions used for the fit.}
	\item{coefficients}{numeric matrix of estimated coefficients from the glm fits. Coefficients are on the natural log scale. Matrix has dimensions \code{nrow(y)} by \code{ncol(design)}.}
	\item{unshrunk.coefficients}{numeric matrix of estimated coefficients from the glm fits when no log-fold-changes shrinkage is applied, on the natural log scale, of size \code{nrow(y)} by \code{ncol(design)}. It exists only when \code{prior.count} is not 0.}
	\item{fitted.values}{matrix of fitted values from glm fits, same number of rows and columns as \code{y}.}
	\item{deviance}{numeric vector of deviances, one for each gene.}
	\item{nthreads}{integer giving the number of OpenMP threads used for the fit, as specified by the \code{nthreads} argument. It is stored on the fitted object so that downstream tests reuse the same number of threads.}
}

\details{
Implements generalized linear model (GLM) methods developed by McCarthy et al (2012).
Specifically, \code{glmFit} fits genewise negative binomial GLMs, all with the same design matrix but possibly different dispersions, offsets and weights.
When the design matrix defines a one-way layout, or can be re-parametrized to a one-way layout, the GLMs are fitted very quickly using Fisher scoring (same algorithm as \code{\link{mglmOneGroup}}).
For more complex design matrices, the GLMs are fitted using Fisher scoring with Levenberg-style damping to ensure convergence (same algorithm as \code{\link{mglmLevenberg}}).

Positive \code{prior.count} values cause the returned coefficients to be shrunk in such a way that fold-changes between the treatment conditions are decreased and infinite fold-changes are avoided (Phipson, 2013).
Larger \code{prior.count} values cause more shrinkage.
Coefficient shrinkage does not affect the likelihood ratio tests or p-values.

Arguments \code{offset}, \code{lib.size} and \code{weights} are not accepted when \code{y} is a DGEList, because these quantities should be defined when the DGEList is created.
}

\references{
McCarthy DJ, Chen Y, Smyth GK (2012).
Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation.
\emph{Nucleic Acids Research} 40, 4288-4297.
\doi{10.1093/nar/gks042}.

Phipson B (2013).
Empirical Bayes modelling of expression profiles and their associations.
Ph.D. thesis, Department of Mathematics and Statistics, The University of Melbourne.
\url{http://hdl.handle.net/11343/38162}.

Chen Y, Chen L, Lun ATL, Baldoni PL, Smyth GK (2025).
edgeR v4: powerful differential analysis of sequencing data with expanded functionality and improved support for small counts and larger datasets.
\emph{Nucleic Acids Research} 53(2), gkaf018.
\doi{10.1093/nar/gkaf018}.
}

\author{Davis McCarthy, Gordon Smyth, Yunshun Chen, Aaron Lun, Lizhong Chen}

\examples{
nlibs <- 3
ngenes <- 100
dispersion.true <- 0.1

# Make first gene respond to covariate x
x <- 0:2
design <- model.matrix(~x)
beta.true <- cbind(Beta1=2,Beta2=c(2,rep(0,ngenes-1)))
linear.predictor.true <- beta.true \%*\% t(design)
mu.true <- 2^linear.predictor.true

# Generate count data
y <- rnbinom(ngenes*nlibs,mu=mu.true,size=1/dispersion.true)
y <- matrix(y,ngenes,nlibs)
colnames(y) <- c("x0","x1","x2")
rownames(y) <- paste("gene",1:ngenes,sep=".")
d <- DGEList(y)

# Normalize
d <- normLibSizes(d)

# Fit the NB GLMs
fit <- glmFit(d, design, dispersion=dispersion.true)

# Likelihood ratio tests for trend
results <- glmLRT(fit, coef=2)
topTags(results)
}

\seealso{
\code{\link{glmQLFit}} is a newer alternative to \code{glmFit} for use in the QL pipeline.

Genewise NB GLMs can also be fitted by the non-generic functions \code{\link{mglmOneGroup}} and \code{\link{mglmLevenberg}}.
Previously these functions were called by \code{glmFit}, but \code{glmFit} now calls C code directly.

\code{\link{topTags}} displays results from \code{glmLRT}.
}

\concept{Model fit}
\concept{Differential expression}
