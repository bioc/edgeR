\name{glmFit}
\alias{glmFit}
\alias{glmFit.DGEList}
\alias{glmFit.SummarizedExperiment}
\alias{glmFit.default}

\title{Genewise Negative Binomial Generalized Linear Models}

\description{Fit a negative binomial generalized log-linear model to the read counts for each gene.}

\usage{
\method{glmFit}{default}(y, design=NULL, dispersion=NULL, offset=NULL, lib.size=NULL, weights=NULL,
       prior.count=0.125, start=NULL, \dots)
\method{glmFit}{DGEList}(y, design=NULL, dispersion=NULL, prior.count=0.125, start=NULL, \dots)
\method{glmFit}{SummarizedExperiment}(y, design=NULL, dispersion=NULL, prior.count=0.125, start=NULL, \dots)
}

\arguments{
\item{y}{a matrix of counts or a \code{DGEList} object or a \code{SummarizedExperiment} containing counts. Rows represent genes and columns represent samples.}

\item{design}{numeric matrix giving the design matrix for the genewise linear models.
Must be of full column rank.
Defaults to a single column of ones, equivalent to treating the columns as replicate libraries.}

\item{dispersion}{numeric scalar, vector or matrix of negative binomial dispersions. Can be a common value for all genes, a vector of dispersion values with one for each gene, or a matrix of dispersion values with one for each observation. If \code{NULL} will be extracted from \code{y}, with order of precedence: genewise dispersion, trended dispersions, common dispersion.}

\item{offset}{numeric matrix of same size as \code{y} giving offsets for the log-linear models.  Can be a scalar or a vector of length \code{ncol(y)}, in which case it is expanded out to a matrix.}

\item{weights}{optional numeric matrix giving prior weights for the observations (for each library and gene) to be used in the GLM calculations.}

\item{lib.size}{numeric vector of length \code{ncol(y)} giving library sizes. Only used if \code{offset=NULL}, in which case \code{offset} is set to \code{log(lib.size)}. Defaults to \code{colSums(y)}.}

\item{prior.count}{average prior count to be added to observation to shrink the estimated log-fold-changes towards zero.}

\item{start}{optional numeric matrix of initial estimates for the linear model coefficients.}

\item{\dots}{other arguments are passed to lower level fitting functions.}
}

\value{
An object of class \code{DGEGLM} containing components \code{counts}, \code{samples}, \code{genes} and \code{abundance} from \code{y} plus the following new components:
	\item{design}{design matrix as input.}
	\item{weights}{matrix of weights as input.}
	\item{df.residual}{numeric vector of residual degrees of freedom, one for each gene.}
	\item{offset}{numeric matrix of linear model offsets.}
	\item{dispersion}{vector of dispersions used for the fit.}
	\item{coefficients}{numeric matrix of estimated coefficients from the glm fits, on the natural log scale, of size \code{nrow(y)} by \code{ncol(design)}.}
	\item{unshrunk.coefficients}{numeric matrix of estimated coefficients from the glm fits when no log-fold-changes shrinkage is applied, on the natural log scale, of size \code{nrow(y)} by \code{ncol(design)}. It exists only when \code{prior.count} is not 0.}
	\item{fitted.values}{matrix of fitted values from glm fits, same number of rows and columns as \code{y}.}
	\item{deviance}{numeric vector of deviances, one for each gene.}
}

\details{
Implements generalized linear model (GLM) methods developed by McCarthy et al (2012).
Specifically, \code{glmFit} fits genewise negative binomial GLMs, all with the same design matrix but possibly different dispersions, offsets and weights.
When the design matrix defines a one-way layout, or can be re-parametrized to a one-way layout, the GLMs are fitting very quickly using \code{\link{mglmOneGroup}}.
Otherwise the default fitting method, implemented in \code{\link{mglmLevenberg}}, uses a Fisher scoring algorithm with Levenberg-style damping.

Positive \code{prior.count} values cause the returned coefficients to be shrunk in such a way that fold-changes between the treatment conditions are decreased and infinite fold-changes are avoided (Phipson, 2013).
Larger \code{prior.count} values cause more shrinkage.
Coefficient shrinkage does not affect the likelihood ratio tests or p-values.
}

\references{
McCarthy DJ, Chen Y, Smyth GK (2012).
Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation.
\emph{Nucleic Acids Research} 40, 4288-4297.
\doi{10.1093/nar/gks042}

Phipson B (2013).
Empirical bayes modelling of expression profiles and their associations.
Ph.D. thesis, Department of Mathematics and Statistics, The University of Melbourne.
\url{http://hdl.handle.net/11343/38162}.

Chen Y, Chen L, Lun ATL, Baldoni PL, Smyth GK (2024). edgeR 4.0: powerful differential analysis of sequencing data with expanded functionality and improved support for small counts and larger datasets.
\emph{bioRxiv} 2024.01.21.576131.
\doi{10.1101/2024.01.21.576131}
}

\author{Davis McCarthy, Gordon Smyth, Yunshun Chen, Aaron Lun}

\examples{
nlibs <- 3
ngenes <- 100
dispersion.true <- 0.1

# Make first gene respond to covariate x
x <- 0:2
design <- model.matrix(~x)
beta.true <- cbind(Beta1=2,Beta2=c(2,rep(0,ngenes-1)))
mu.true <- 2^(beta.true \%*\% t(design))

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
Low-level computations are done by \code{\link{mglmOneGroup}} or \code{\link{mglmLevenberg}}.

\code{\link{topTags}} displays results from \code{glmLRT}.
}

\concept{Model fit}
\concept{Differential expression}
