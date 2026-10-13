\name{filterByExpr}
\alias{filterByExpr}
\alias{filterByExpr.DGEList}
\alias{filterByExpr.SummarizedExperiment}
\alias{filterByExpr.default}

\title{Filter Genes by Expression Level}

\description{Determine which genes have sufficiently large counts to be retained in a statistical analysis.}

\usage{
\method{filterByExpr}{default}(
#    Data matirx
     y,
    
#    Minimum sample size
     design = NULL,
     group = NULL,
     min.sample.size = NULL,

#    Library sizes
     lib.size = NULL,

#	 Minimum counts per sample and in total
     min.count = 10,
     min.total.count = 15,

#    Limit the sample size for large group sizes
     large.n = 10,
     min.prop = 0.7,
     max.n.samples = 40,

#    Preset options
     preset = NULL,
     
#    Any other arguments    
     \dots)

\method{filterByExpr}{DGEList}(y, design = NULL, group = NULL, lib.size = NULL, \dots)
\method{filterByExpr}{SummarizedExperiment}(y, design = NULL, group = NULL, lib.size = NULL, \dots)
}

\arguments{ 
\item{y}{matrix of counts, or a \code{DGEList} object, or a \code{SummarizedExperiment} object.}
\item{design}{
  design matrix.
  Defaults to \code{y$design} if \code{y} is a DGEList.
  Used to compute the minimum group sample size (\code{MinSampleSize}).
  Ignored if \code{group} or \code{min.sample.size} are not \code{NULL}.
}
\item{group}{
  vector or factor giving group membership for the samples.
  Provides an alternative way to specify the design matrix when there is only one explanatory factor.
  Defaults to \code{y$samples$group} if \code{y} is a DGEList and \code{design} and \code{group} are both \code{NULL}.
  Ignored if \code{min.sample.size} are not \code{NULL}.
}
\item{min.sample.size}{
  minimum group sample size.
  Provides a way to specify the minimum group sample size directly, instead of allowing it to be computed from the \code{design} or \code{group} arguments.
  If not \code{NULL}, then takes precedence over \code{design} or \code{group} for specifying \code{MinSampleSize}.
}
\item{lib.size}{
  vector of library sizes.
  Defaults to \code{colSums(y)} if \code{y} is a matrix or to \code{normLibSizes(y)} if \code{y} is a DGEList.
  Used to compute the CPM cutoff value.
}
\item{min.count}{
  minimum read count required for at least some samples.
  Used to compute the CPM cutoff value.
}
\item{min.total.count}{
  minimum total count required across all samples.
}
\item{large.n}{
  number of samples that would be considered a \dQuote{large} group.
  Used to modify \code{MinSampleSize} for large datasets.
}
\item{min.prop}{
  proportion of samples in large groups for which a gene needs to be expressed to be biologically interesting.
  Used to modify \code{MinSampleSize} for large datasets.
  See Details below for the exact formula.
}
\item{max.n.samples}{
  upper bound on the number of samples in which a gene needs to be expressed.
}
\item{preset}{
  choose a preset combination of tuning arguments.
  Possibilities are \code{"stringent"}, \code{"moderate"}, or \code{"lenient"}, going successively from heavy to light filtering.
  The default settings match the \code{"stringent"} option.
  If this argument is used then it overrides any settings made for individual tuning arguments.
}
\item{\dots}{
  other arguments are not currently used.
}
}

\details{
This function implements the filtering strategy that was described informally by Chen et al (2016).
Roughly speaking, the strategy keeps genes that have at least \code{min.count} reads in a worthwhile number of samples.

More precisely, the filtering keeps genes that have \code{CPM >= CPM.cutoff} in \code{MinSampleSize} samples,
where \code{CPM.cutoff = min.count / median(lib.size) * 1e6} and \code{MinSampleSize} is the smallest group sample size or, more generally, the minimum inverse leverage computed from the design matrix.

If all the group sample sizes are large, then \code{MinSampleSize} is reduced somewhat.
If \code{MinSampleSize > large.n}, then genes are kept if \code{CPM >= CPM.cutoff} in \code{k} samples where
\code{k = min( large.n + (MinSampleSize - large.n) * min.prop, max.n.samples)}.
This rule requires that genes are expressed in at least \code{min.prop * MinSampleSize} samples when \code{MinSampleSize} is large, with \code{max.n.samples} being the upper bound on the number of samples that a gene much be expressed in.

In addition, each kept gene is required to have at least \code{min.total.count} reads across all the samples.

In edgeR v5, the new arguments \code{min.sample.size}, \code{max.n.samples}, and \code{preset} have been added.
\code{min.sample.size} provides a way to specify the minimum effective group size for the experiment without inputing \code{design} or \code{group}.
\code{max.n.samples} puts an upper limit on the number of samples in which a gene must be expressed to be kept in the analysis, a consideration that becomes relevant for large-scale experiments with many samples in each group.

The \code{"preset"} argument allows users to choose a preset combination of tuning parameters instead of setting arguments individually.
Choose \code{"stringent"} to keep only genes with enough reads across enough samples to have a good chance of being called as DE if they are genuinely DE.
This choice corresponds to the historical default settings of \code{filterByExpr}.
Choose \code{"lenient"} to keep the maximum number of genes in the analysis, including all genes that could possibly be called as DE if the true differential differences are large and the dispersion levels are not too high.
\code{"moderate"} provides an intermediate setting between the two.
}

\note{
Note that \code{filterByExpr} implements independent filtering in the sense that the filtering would be exactly the same even if the sample group labels or the rows of the design matrix were permuted. 
The function uses the \code{design} and \code{group} arguments only to identify the smallest group size (\code{MinSampleSize}). 
Non-independent filtering that used information about which sample belongs to each group, for example keeping genes expressed in a certain number of samples in at least one group, would be invalid because it would potentially cause false discoveries in subsequent differential expression analyses. 
}

\value{
Logical vector of length \code{nrow(y)} indicating which rows of \code{y} to keep in the analysis.
}

\author{Gordon Smyth}

\references{
Chen Y, Lun ATL, Smyth GK (2016).
From reads to genes to pathways: differential expression analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline.
\emph{F1000Research} 5, 1438.
\doi{10.12688/f1000research.8987.2}.

Law CW, Alhamdoosh M, Su S, Dong X, Tian L, Smyth GK, Ritchie ME (2016).
RNA-seq analysis is easy as 1-2-3 with limma, Glimma and edgeR.
\emph{F1000Research} 5, 1408.
\doi{https://doi.org/10.12688/f1000research.9005.3}.
}

\examples{
condition <- factor(c(1,1,1,2,2,2))
levels(condition) <- c("A","B")
counts <- matrix(0,4,6)
row.names(counts) <- paste0("Gene",1:4)
colnames(counts) <- paste0(condition,c(1,2,3,1,2,3))
counts[1,c(4,5,6)] <- 7
counts[2,c(4,5)] <- 8
counts[3,c(2,3,5)] <- 8
counts[4,c(4,5,6)] <- 8
y <- DGEList(counts=counts,lib.size=rep(1e6,6))
keep <- filterByExpr(y, group=condition, min.count=8)
data.frame(counts,Keep=keep)
}

\concept{Preprocessing}
