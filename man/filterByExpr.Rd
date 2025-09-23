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
\item{design}{
  design matrix.
  Ignored if \code{group} is not \code{NULL}.
  Defaults to \code{y$design} if \code{y} is a DGEList.
  Used to compute the minimum group sample size (\code{MinSampleSize}).
}
\item{group}{
  vector or factor giving group membership.
  Provides an alternative way to specify the design matrix when there is only one explanatory factor (oneway layout).
  Defaults to \code{y$samples$group} if \code{y} is a DGEList and \code{design} and \code{group} are both \code{NULL}.
}
\item{lib.size}{
  library size.
  Defaults to \code{colSums(y)} if \code{y} is a matrix or to \code{normLibSizes(y)} if \code{y} is a DGEList.
  Used to compute the CPM cutoff value.
}
\item{min.count}{
  numeric.
  Minimum count required for at least some samples.
  Used to compute the CPM cutoff value.
}
\item{min.total.count}{numeric. Minimum total count required across all samples.}
\item{large.n}{
  integer.
  Number of samples that would be considered a \dQuote{large} group.
  Used to modify \code{MinSampleSize} for large datasets.
}
\item{min.prop}{
  numeric.
  Used to modify \code{MinSampleSize} for large datasets.
  See Details below for the exact formula.
  }
\item{\dots}{other arguments are not currently used.}
}

\details{
This function implements the filtering strategy that was described informally by Chen et al (2016).
Roughly speaking, the strategy keeps genes that have at least \code{min.count} reads in a worthwhile number samples.

More precisely, the filtering keeps genes that have \code{CPM >= CPM.cutoff} in \code{MinSampleSize} samples,
where \code{CPM.cutoff = min.count / median(lib.size) * 1e6} and \code{MinSampleSize} is the smallest group sample size or, more generally, the minimum inverse leverage computed from the design matrix.

If all the group samples sizes are large, then \code{MinSampleSize} is reduced somewhat.
If \code{MinSampleSize > large.n}, then genes are kept if \code{CPM >= CPM.cutoff} in \code{k} samples where
\code{k = large.n + (MinSampleSize - large.n) * min.prop}.
This rule requires that genes are expressed in at least \code{min.prop * MinSampleSize} samples, even when \code{MinSampleSize} is large.

In addition, each kept gene is required to have at least \code{min.total.count} reads across all the samples.
}

\note{
Note that \code{filterByExpr} implements independent filtering in the sense that the filtering would be exactly the same even if the sample group labels or the rows of the design matrix were permuted. 
The function uses the \code{design} and \code{group} arguments only to identify the smallest group size (\code{MinSampleSize}). 
Non-independent filtering that used information about which sample belongs to each group, for example keeping genes expressed in a certain nmber of samples in at least one group, would be invalid because it would potentially cause false discoveries in subsequent differential expression analyses. 
}

\value{
Logical vector of length \code{nrow(y)} indicating which rows of \code{y} to keep in the analysis.
}

\author{Gordon Smyth}

\references{
Chen Y, Lun ATL, and Smyth, GK (2016).
From reads to genes to pathways: differential expression analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline.
\emph{F1000Research} 5, 1438.
\doi{10.12688/f1000research.8987.2}
}

\examples{\dontrun{
keep <- filterByExpr(y, design)
y <- y[keep,]
}}
