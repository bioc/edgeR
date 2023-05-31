\name{filterByExpr}
\alias{filterByExpr}
\alias{filterByExpr.DGEList}
\alias{filterByExpr.SummarizedExperiment}
\alias{filterByExpr.default}

\title{Filter Genes By Expression Level}

\description{Determine which genes have sufficiently large counts to be retained in a statistical analysis.}

\usage{
\method{filterByExpr}{DGEList}(y, design = NULL, group = NULL, lib.size = NULL, \dots)
\method{filterByExpr}{SummarizedExperiment}(y, design = NULL, group = NULL, lib.size = NULL, \dots)
\method{filterByExpr}{default}(y, design = NULL, group = NULL, lib.size = NULL,
             min.count = 10, min.total.count = 15, large.n = 10, min.prop = 0.7, \dots)
}

\arguments{ 
\item{y}{matrix of counts, or a \code{DGEList} object, or a \code{SummarizedExperiment} object.}
\item{design}{design matrix. Ignored if \code{group} is not \code{NULL}.}
\item{group}{vector or factor giving group membership for a oneway layout, if appropriate.}
\item{lib.size}{library size, defaults to \code{colSums(y)}.}
\item{min.count}{numeric. Minimum count required for at least some samples.}
\item{min.total.count}{numeric. Minimum total count required.}
\item{large.n}{integer. Number of samples per group that is considered to be \dQuote{large}.}
\item{min.prop}{numeric. In large sample situations, the minimum proportion of samples in a group that a gene needs to be expressed in. See Details below for the exact formula.}
\item{\dots}{any other arguments.
For the \code{DGEList} and \code{SummarizedExperiment} methods, other arguments will be passed to the default method.
For the default method, other arguments are not currently used.}
}

\details{
This function implements the filtering strategy that was described informally by Chen et al (2016).
Roughly speaking, the strategy keeps genes that have at least \code{min.count} reads in a worthwhile number samples.

More precisely, the filtering keeps genes that have \code{CPM >= CPM.cutoff} in \code{MinSampleSize} samples,
where \code{CPM.cutoff = min.count/median(lib.size)*1e6} and \code{MinSampleSize} is the smallest group sample size or, more generally, the minimum inverse leverage computed from the design matrix.

If all the group samples sizes are large, then the above filtering rule is relaxed slightly.
If \code{MinSampleSize > large.n}, then genes are kept if \code{CPM >= CPM.cutoff} in \code{k} samples where
\code{k = large.n + (MinSampleSize - large.n) * min.prop}.
This rule requires that genes are expressed in at least \code{min.prop * MinSampleSize} samples, even when \code{MinSampleSize} is large.

In addition, each kept gene is required to have at least \code{min.total.count} reads across all the samples.
}

\value{
Logical vector of length \code{nrow(y)} indicating which rows of \code{y} to keep in the analysis.
}

\author{Gordon Smyth}

\references{
Chen Y, Lun ATL, and Smyth, GK (2016).
From reads to genes to pathways: differential expression analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline.
\emph{F1000Research} 5, 1438.
\url{https://f1000research.com/articles/5-1438}
}

\examples{\dontrun{
keep <- filterByExpr(y, design)
y <- y[keep,]
}}
