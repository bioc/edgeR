\name{glmfit}
\alias{glmFit}
\alias{glmLRT}

\title{Fit negative binomial generalized linear model for each transcript}

\description{Fit a negative binomial generalized linear model for each transcript (tag) with the unadjusted counts provided, a value for the dispersion parameter and, optionally, offsets and weights for different libraries or transcripts.}


\usage{glmFit(y, design, dispersion, offset=0, weights=NULL, lib.size=NULL) 
glmLRT(y, glmfit, coef=ncol(glmfit$design))}

\arguments{ 

\item{y}{an object that contains the raw counts for each library (the measure of expression level); alternatively, a matrix of counts, or a \code{DGEList} object with (at least) elements \code{counts} (table of unadjusted counts) and \code{samples} (data frame containing information about experimental group, library size and normalization factor for the library size)}

\item{design}{numeric matrix giving the design matrix for the GLM that is to be fit.}

\item{dispersion}{numeric scalar or vector providing the value for the dispersion parameter that is used in fitting the GLM for each transcript. Can be a common value for all tags, or a vector of values can provide a unique dispersion value for each tag.}

\item{offset}{numeric scalar, vector or matrix giving the offset (in addition to the log of the effective library size) that is to be included in the NB GLM for the transcripts. If a scalar, then this value will be used as an offset for all transcripts and libraries. If a vector, it should be have length equal to the number of libraries, and the same vector of offsets will be used for each transcript. If a matrix, then each library for each transcript can have a unique offset, if desired.}

\item{weights}{optional numeric matrix giving the matrix of weights for the observations (for each library and transcript) to be used in the GLM calculations.}

\item{lib.size}{optional vector providing the (effective) library size for each library (must have length equal to the number of columns, or libraries, in the matrix of counts). If \code{NULL}, then a default is used. If \code{y} is a \code{DGEList} object then the default for \code{lib.size} is the product of the library sizes and the normalization factors (in the \code{samples} slot of the object). If \code{y} is a simple matrix of counts, then the default for \code{lib.size} is the vector of column sums of \code{y}.}

\item{glmfit}{a \code{DGEGLM} object, the output from \code{glmFit}.}

\item{coef}{scalar or vector indicating the column(s) of \code{design} that are to be dropped when creating the null model for the Likelihood Ratio (LR) Test. The \code{glmLRT} fits the null model and then conducts an LR test of the model fit provided in \code{glmfit} against the null model defined by the choice of \code{coef}.}
}

\value{\code{glmFit} produces an object of class \code{DGEGLM} with the following components:
	\item{coefficients}{matrix of estimated coefficients from the NB model}
	\item{df.residual}{vector giving the residual degrees of freedom for each tag. In theory it can be different for different tags (if there are missing values), but in practice these will usually be identical for each tag.}
	\item{deviance}{vector giving the deviance from the NB model fit for each tag.}
	\item{design}{design matrix used in the NB model fit for each tag.}
	\item{offset}{scalar, vector or matrix giving the offset to use in the NB model for each tag.}
	\item{samples}{data frame providing information about the samples (libraries) in the experiment; taken from the object \code{y}.}
	\item{genes}{vector or data frame providing gene information for each tag; taken from the object \code{y}.}
        \item{dispersion}{scalar or vector giving the the value of the dispersion parameter used in each tag's NB model fit.}
        \item{lib.size}{vector of library sizes used in the model fit.}
        \item{weights}{matrix of final weights used in the NB model fits for each tag.}
        \item{fitted.values}{matrix of fitted values from the NB model for each tag.}
        \item{abundance}{vector of gene/tag abundances (expression level), on the log2 scale, computed from the mean count for each gene/tag after scaling count by normalized library size.}

\code{glmLRT} produces an object of class \code{DGELRT} with the following components:
              \item{table}{data frame (table) containing the abundance of each tag (log-concentration, \code{logConc}), the log-fold change of expression between conditions/contrasts being tested (\code{logFC}), the likelihood ratio statistic (\code{LR.statistic}) and the p-value from the LR test (\code{p.value}), for each tag in the dataset.}
              \item{coefficients}{matrix of coefficients for the full model defined by the \code{design} matrix (i.e. for the full model).}
              \item{dispersion.used}{scalar or vector of the dispersion value(s) used in the GLM fits and LR test.}
The \code{DGELRT} object also contains all the elements of \code{y} except for the table of counts (raw data) and the table of pseudo-counts (if applicable).

}

\details{
Given a fixed value for the dispersion parameter, a negative binomial model can be fitted to the counts for each tag/transcript in a dataset. The function \code{glmFit} calls the in-built function \code{glm.fit} to fit the NB GLM for each tag. Once we have a fit for a given design matrix, \code{glmLRT} can be run with a given coefficient or contrast specified and evidence for differential expression assessed using a likelihood ratio test. Tags can be ranked in order of evidence for differential expression, based on the p-value computed for each tag.
}


\author{Davis McCarthy and Gordon Smyth}
\examples{
nlibs <- 3
ntags <- 100
dispersion.true <- 0.1

# Make first transcript respond to covariate x
x <- 0:2
design <- model.matrix(~x)
beta.true <- cbind(Beta1=2,Beta2=c(2,rep(0,ntags-1)))
mu.true <- 2^(beta.true \%*\% t(design))

# Generate count data
y <- rnbinom(ntags*nlibs,mu=mu.true,size=1/dispersion.true)
y <- matrix(y,ntags,nlibs)
colnames(y) <- c("x0","x1","x2")
rownames(y) <- paste("Gene",1:ntags,sep="")
d <- DGEList(y)

# Normalize
d <- calcNormFactors(d)

# Fit the NB GLMs
fit <- glmFit(d, design, dispersion=dispersion.true)

## Likelihood ratio tests for trend
results <- glmLRT(d, fit, coef=2)
topTags(results)
}

\seealso{
\code{\link{estimateCRDisp}} for estimating the negative binomial dispersion.

\code{\link{topTags}} for displaying results from \code{glmLRT}.
}

\keyword{algebra}