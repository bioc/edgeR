\name{glmLRT}
\alias{glmLRT}

\title{Genewise Likelihood Ratio Tests}

\description{Given genewise generalized linear model fits, conduct likelihood ratio tests for a given coefficient or coefficient contrast.}

\usage{
glmLRT(glmfit, coef=ncol(glmfit$design), contrast=NULL)
}

\arguments{
\item{glmfit}{a \code{DGEGLM} object, usually output from \code{glmFit}.}

\item{coef}{integer or character vector indicating which coefficients of the linear model are to be tested equal to zero.  Values must be columns or column names of \code{design}. Defaults to the last coefficient.  Ignored if \code{contrast} is specified.}

\item{contrast}{numeric vector or matrix specifying one or more contrasts of the linear model coefficients to be tested equal to zero. If a matrix, then rows correspond to the columns of \code{design} and columns are contrasts. If non-NULL, then takes precedence over \code{coef}.}
}

\value{
An object of class \code{DGELRT} with the same components as for \code{\link{glmFit}} plus the following:
	\item{table}{data frame with the same rows as \code{y} containing the log2-fold-changes, likelhood ratio statistics and p-values, ready to be displayed by \code{topTags}.}
	\item{comparison}{character string describing the coefficient or the contrast being tested.}

The data frame \code{table} contains the following columns:
	\item{logFC}{log2-fold change of expression between conditions being tested.}
	\item{logCPM}{average log2-counts per million, the average taken over all libraries in \code{y}.}
	\item{LR}{likelihood ratio statistics.}
	\item{PValue}{p-values.}
}

\details{
Using genewise GLMs from \code{glmFit}, \code{glmLRT} conducts likelihood ratio tests for one or more coefficients in the linear model.
If \code{contrast} is NULL, the null hypothesis is that all the coefficients specified by \code{coef} are equal to zero.
If \code{contrast} is non-NULL, then the null hypothesis is that the specified contrasts of the coefficients are equal to zero.
For example, \code{contrast = c(0,1,-1)}, assuming there are three coefficients, would test the hypothesis that the second and third coefficients are equal.
If \code{contrast} is a matrix, then each column is a contrast and the null hypothesis is that all the contrasts are equal to zero.
}

\references{
McCarthy DJ, Chen Y, Smyth GK (2012). Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation.
\emph{Nucleic Acids Research} 40, 4288-4297.
\doi{10.1093/nar/gks042}
}

\author{Gordon Smyth, Davis McCarthy, Yunshun Chen}

\examples{
# See glmFit for an example that includes glmLRT
}

\seealso{
\code{\link{glmFit}} fits the genewise GLMs that are input to \code{glmLRT}.

\code{\link{topTags}} displays results from \code{glmLRT}.
}

\concept{Model fit}
\concept{Differential expression}
