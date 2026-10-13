\name{splitGencodeTxNames}
\alias{splitGencodeTxNames}
\alias{collapseGencodeTx2Gene}

\title{Split GENCODE Transcript Names and Summarize By Gene}

\description{Split delimited transcript names provided by GENCODE FASTA files into separate columns.}

\usage{
splitGencodeTxNames(x, remove.version = TRUE)
collapseGencodeTx2Gene(txann)
}

\arguments{
\item{x}{character vector of GENCODE transcript names.}
\item{remove.version}{logical, whether to remove version numbers from transcript and gene IDs.}
\item{txann}{eight-column character matrix produced by \code{splitGencodeTxNames}.}
}

\details{
GENCODE transcript names consist of different annotations pasted together with \code{|} delimiters.
\code{splitGencodeTxNames} splits the annotations into an eight-column character matrix.
\code{collapseGencodeTx2Gene} collapses the output of \code{splitGencodeTxNames} to a gene-level data.frame.
}

\value{
\code{splitGencodeTxNames} produces a character matrix with the following columns:
  \item{EnsT}{Ensembl transcript ID, without version number suffix.}
  \item{EnsG}{Ensembl gene ID, without version number suffix.}
  \item{HavG}{Havana gene ID.}
  \item{HavT}{Havana transcript ID.}
  \item{TxName}{Transcript name.}
  \item{GeneName}{Gene name.}
  \item{AnnLength}{Transcript length.}
  \item{Type}{Biotype of transcript.}

\code{collapseGencodeTx2Gene} produces a data.frame with Ensembl gene IDs as row.names and the following columns:
  \item{GeneName}{Gene name.}
  \item{TxTypes}{Unique transcript biotypes for that gene (pipe-separated).}
  \item{NPCTx}{Number of protein-coding transcripts for that gene.}
}

\author{Gordon Smyth}

\seealso{
\code{\link{catchSalmon}}, \code{\link{catchSalmonGene}}.
}

\concept{Transcript expression}
\concept{Reading data files}
