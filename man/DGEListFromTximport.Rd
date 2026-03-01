\name{DGEListFromTximport}
\alias{DGEListFromTximport}

\title{
Create a DGEList Object From tximport Output
}

\description{
Assembles a \code{DGEList} object from the list object created by \code{tximport::tximport}.
}

\usage{
DGEListFromTximport(txi, samples = NULL, group = NULL, genes = NULL,
    remove.zeros = FALSE, norm.method = "none")
}

\arguments{
  \item{txi}{list output from \code{tximport}.}
  \item{samples}{data.frame containing sample information, with a row for each sample. This data.frame will be appended to the \code{samples} component of the DGEList object.}
  \item{group}{vector or factor giving the experimental group or treatment condition for each sample. Defaults to a single group.}
  \item{genes}{data.frame containing gene annotation.}
  \item{remove.zeros}{logical, whether to remove rows that have 0 total count.}
  \item{norm.method}{library size normalization method, if desired. See \code{normLibSizes} for options.}
}

\details{
This function imports the output \code{txi} from the \code{tximport::tximport} function to the edgeR package.
In most use cases, \code{tximport}, will have been run with \code{txOut=FALSE} to aggregate counts at the level of genes, and the transcript length matrix \code{txi$length} will contain the average transcript length for each gene in each sample.
Alternatively, \code{txi} will contain results at the original transcript level, if \code{tximport} was run with \code{txOut=TRUE}, in which case \code{txi$length} will contain effective transcript lengths for the individual transcripts in each sample. 

The count matrix \code{txi$counts} becomes the count matrix for the DGEList object.
If \code{txi$countsFromAbundance=="no"}, then a prior offset matrix is formed from the average transcript lengths found in \code{txi$length}.
This implements the "original counts and offset" method described by the tximport documentation.
If \code{txi$countsFromAbundance} has some other value, then the counts have already be bias corrected and no offset matrix will be formed.
This implements the "bias corrected counts without an offset" method described by the tximport documentation.

If \code{txi$infReps} exists and contains technical replicates (boostrap or Gibbs samples) for each column, then an read-assignment-ambiguity overdispersion is estimated for each row of data.

Library size normalization will be done if \code{norm.method} specifies a method offered by \code{normLibSizes}, but performing library size normalization as part of the import is not necessary.
We generally recommend leaving \code{norm.method="none"} and running \code{normLibSizes} on the DGEList later in the analysis pipeline, typically after filtering.
}

\value{
A \code{\link[edgeR:DGEList-class]{DGEList}} object.
If \code{txi$counts} contains raw counts, then the a prior offset matrix \code{offset.prior} will be included to account for variation in average transcript length between samples.
The \code{genes} component of the object will contain the following columns:
  \item{TxLengthAve}{average transcript length, over all samples and all transcripts.}
  \item{TxLengthRange}{ratio of largest to smallest transcript length for different samples. If counts are at the gene level, then the transcript lengths are averages over transcripts for each sample.}
  \item{Overdispersion}{quasi-dispersion arising from read assignment ambiguity. Equal to 1 for rows with no ambiguity.}
}

\references{
Baldoni PL, Chen Y, Hediyeh-zadeh S, Liao Y, Dong X, Ritchie ME, Shi W, Smyth GK (2024a).
Dividing out quantification uncertainty allows efficient assessment of differential transcript expression with edgeR.
\emph{Nucleic Acids Research} 52(3), e13.
\doi{10.1093/nar/gkad1167}.

Baldoni PL, Chen L, Li M, Chen Y, Smyth GK (2025).
Dividing out quantification uncertainty enables assessment of differential transcript usage with limma and edgeR.
\emph{Nucleic Acids Research} 53(22), gkaf1305.
\doi{10.1093/nar/gkaf1305}
}

\seealso{
\code{\link{catchSalmon}}, \code{\link[edgeR:DGEList-class]{DGEList-class}}.

The input object \code{txi} is created by \code{tximport} in the tximport package.
}

\concept{edgeR classes}
\concept{Reading data files}
