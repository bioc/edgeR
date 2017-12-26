\name{filterByExpr}
\alias{filterByExpr}
\alias{filterByExpr.DGEList}
\alias{filterByExpr.default}

\title{Filter Genes By Expression Level}

\description{Determine which genes have sufficiently large counts to be retained in a statistal analysis.}

\usage{
\method{filterByExpr}{DGEList}(y, design = NULL, group = NULL, lib.size = NULL, \dots)
\method{filterByExpr}{default}(y, design = NULL, group = NULL, lib.size = NULL,
             min.count = 10, min.total.count = 15, \dots)
}

\arguments{ 
\item{y}{matrix of counts or a \code{DGEList} object.}
\item{design}{design matrix. Ignored if \code{group} is not \code{NULL}.}
\item{group}{vector or factor giving group membership for a oneway layout, if appropriate.}
\item{lib.size}{library size, defaults to \code{colSums(y)}.}
\item{min.count}{numeric. Minimum count required for at least some samples.}
\item{min.total.count}{numeric. Minimum total count required.}
\item{\dots}{any other arguments.
For the \code{DGEList} methods, other arguments will be passed to the default method.
For the default method, other arguments are not currently used.}
}

\value{
Logical vector of length \code{nrow(y)} indicating which rows of \code{y} to keep in the analysis.
}

\author{Gordon Smyth}
