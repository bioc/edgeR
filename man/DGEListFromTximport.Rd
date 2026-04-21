\name{DGEListFromTximport}
\alias{DGEListFromTximport}
\alias{DGEListFromTximeta}

\title{
Create a DGEList Object From tximport or tximeta Output
}

\description{
Assembles a \code{DGEList} object from the list object created by either \code{tximport::tximport} or \code{tximeta::tximeta}.
}

\usage{
DGEListFromTximport(txi, samples = NULL, group = NULL, genes = NULL,
    remove.zeros = FALSE, divide = FALSE)
DGEListFromTximeta(txm, samples = NULL, group = NULL,
    remove.zeros = FALSE, divide = FALSE)
}

\arguments{
  \item{txi}{list output from \code{tximport()}.}
  \item{txm}{SummarizedExperiment data object from \code{tximeta()}.}
  \item{samples}{data.frame containing sample information, with a row for each sample. This data.frame will be appended to the \code{samples} component of the DGEList object.}
  \item{group}{vector or factor giving the experimental group or treatment condition for each sample. Defaults to a single group.}
  \item{genes}{data.frame containing gene annotation.}
  \item{remove.zeros}{logical, whether to remove rows that have 0 total count.}
  \item{divide}{logical, whether to divide the counts by the RAA overdispersions.}
}

\details{
The functions \code{DGEListFromTximport} and \code{DGEListFromTximeta} import output from the \code{tximport::tximport} and \code{tximeta::tximeta} functions, respectively, to the edgeR package.
In many cases, \code{tximport} and \code{tximeta} will have been run with \code{txOut=FALSE} to aggregate counts at the gene level.
In that case, the length matrix (\code{txi$length} for \code{tximport} and \code{assay(txi,"length")} for \code{tximeta}) will contain the average transcript length for each gene in each sample.
Alternatively, \code{txi} or \code{txm} will contain results at the original transcript level, if \code{tximport} and \code{tximeta} were run with \code{txOut=TRUE}, in which case the length matrix will contain effective transcript lengths for the individual transcripts in each sample. 

The count matrix, \code{txi$counts} for \code{tximport} and \code{assay(txi,"counts")} for \code{tximeta}, becomes the count matrix for the DGEList object.
If either \code{tximport} or \code{tximeta} were run with the option \code{countsFromAbundance="no"}, then a prior offset matrix is formed from the average transcript lengths found in the respective transcript length matrix.
This implements the "original counts and offset" method described by the tximport documentation.
If \code{countsFromAbundance} has been set to any other value, then the counts have already be bias corrected and no offset matrix will be formed.
This implements the "bias corrected counts without an offset" method described by the tximport documentation.

If inferential replicates (bootstrap or Gibbs resampling) exist in the \code{txi} or \code{txm} data objects then an read-assignment-ambiguity (RAA) overdispersion is estimated for each row of data, using the approach described by Baldoni et al (2024a).
If \code{divide=TRUE}, then divided counts will be formed by dividing the original counts by the RAA overdispersions.

Note that transcript level counts from various quantification tools can also be read directly by the native edgeR functions
\code{\link{catchSalmon}}, \code{\link{catchKallisto}}, \code{\link{catchRSEM}}, or \code{\link{catchOarfish}}, without the need for \code{tximport} as an intermediary, and these will generally give the same results as \code{DGEListFromTximport} at the transcript level.
In the case of RSEM, using \code{catchRSEM} is preferable, because RSEM does not output inference replicates but \code{catchRSEM} can nevertheless estimate the assignment over-dispersion from the inferential replicate means and variances.
At the gene level, \code{tximport} and \code{DGEListFromTximport} implement transcript length normalization, which the direct edgeR functions do not.
}

\value{
A \code{\link[edgeR:DGEList-class]{DGEList}} object.

If \code{txi} or \code{txm} contains raw counts, then the a prior offset matrix \code{offset.prior} will be included to account for variation in average transcript length between samples.

The \code{genes} component of the object will contain the following columns:
  \item{AveLength}{average transcript or gene length, over all samples and all transcripts.}
  \item{Max2MinLength}{ratio of largest to smallest length for different samples. If counts are at the gene level, then the lengths are averages over transcripts for each sample, as computed by tximport.}
  \item{Overdispersion}{quasi-dispersion arising from read assignment ambiguity. Equal to 1 for rows with no ambiguity. Larger values indicate more ambiguity.}

As well as standard DGEList components, the output components \code{divided.counts} indicates whether the counts have been divided and the component \code{txi.counts} indicates whether the tximported counts were direct from the quantifier or whether they were scaled by TPMs.
}

\references{
Baldoni PL, Chen Y, Hediyeh-zadeh S, Liao Y, Dong X, Ritchie ME, Shi W, Smyth GK (2024a).
Dividing out quantification uncertainty allows efficient assessment of differential transcript expression with edgeR.
\emph{Nucleic Acids Research} 52(3), e13.
\doi{10.1093/nar/gkad1167}.

Baldoni PL, Chen L, Smyth GK (2024b).
Faster and more accurate assessment of differential transcript expression with Gibbs sampling and edgeR v4.
\emph{NAR Genomics and Bioinformatics} 6(4), lqae151.
\doi{10.1093/nargab/lqae151}.

Baldoni PL, Chen L, Li M, Chen Y, Smyth GK (2025).
Dividing out quantification uncertainty enables assessment of differential transcript usage with limma and edgeR.
\emph{Nucleic Acids Research} 53(22), gkaf1305.
\doi{10.1093/nar/gkaf1305}
}

\seealso{
\code{\link{catchSalmon}}, \code{\link{catchKallisto}}, \code{\link{catchRSEM}}, \code{\link{catchOarfish}}, \code{\link[edgeR:DGEList-class]{DGEList-class}}.

The input objects \code{txi} or \code{txm} are created by either \code{tximport} or \code{tximeta} in the tximport and tximeta packages, respectively.
}

\concept{edgeR classes}
\concept{Reading data files}
