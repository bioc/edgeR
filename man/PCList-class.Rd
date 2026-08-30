\name{PCList-class}
\alias{PCList-class}
\docType{class}
\title{Paired Counts Data Object - Class}

\description{
A list-based S4 class for storing paired counts and associated information from sequencing technologies.
}

\section{List Components}{
For objects of this class, rows correspond to genomic features and columns to samples.
The genomic features can be genes, genomic loci, or any feature for which paired read counts are meaningful.
Objects of this class contain the following essential list components:
\describe{
  \item{\code{counts}:}{ numeric matrix of read counts, one row for each feature and one column for each sample.}
  \item{\code{counts2}:}{ numeric matrix of paired read counts, one row for each feature and one column for each sample.}
  \item{\code{samples}:}{ data.frame with a row for each sample and column \code{group} containing the group labels.
  Other columns can be optionally added to give more detailed sample information.}
}
Optional components include:
\describe{
  \item{\code{genes}:}{ data.frame giving annotation information for each feature.  Same number of rows as \code{counts}.}
  \item{\code{AveLogCPM}:}{ numeric vector giving average log2 total counts per million for each feature.}
}
}

\seealso{
\code{\link{PCList}} constructs PCList objects.
Other classes defined in edgeR are \code{\link{DGEList-class}}, \code{\link{DGEExact-class}}, \code{\link{DGEGLM-class}}, \code{\link{DGELRT-class}}, \code{\link{TopTags-class}}
}

\section{Methods}{
This class inherits directly from class \code{list}, so \code{PCList} objects can be manipulated as if they were ordinary lists.
However they can also be treated as if they were matrices for the purposes of subsetting.

The dimensions, row names and column names of a \code{PCList} object are defined by those of \code{counts}, see \code{\link{dim.DGEList}} or \code{\link{dimnames.DGEList}}.

\code{PCList} objects can be subsetted, see \code{\link[edgeR:subsetting]{subsetting}}.

\code{PCList} objects inherit a compact \code{show} method from the \code{LargeDataObject} class, so that printing produces a compact summary of their contents.
}

\author{edgeR team. First created by Lizhong Chen.}

\concept{edgeR classes}

