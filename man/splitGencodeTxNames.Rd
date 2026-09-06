\name{splitGencodeTxNames}
\alias{splitGencodeTxNames}

\title{Split Gencode Transcript Names}

\description{Split delimited transcript names provided by Gencode FASTA files into separate columns.}

\usage{
splitGencodeTxNames(x, remove.version.numbers = TRUE)
}

\arguments{
\item{x}{character vector of Gencode transcript names.}
\item{remove.version.numbers}{logical, whether to remove version numbers from transcript and gene IDs.}
}

\details{
Gencode transcript names consist of different annotations pasted together with \code{|} delimiters.
Here we split the annotations into separate columns of a matrix.
}

\value{
A character matrix with the following columns:
  \item{EnsT}{Ensembl transcript ID, without version number suffix.}
  \item{EnsG}{Ensembl gene ID, without version number suffix.}
  \item{HavG}{Havana gene ID.}
  \item{HavT}{Havana transcript ID.}
  \item{TxName}{Transcript name.}
  \item{GeneName}{Gene name.}
  \item{AnnLength}{Transcript length.}
  \item{Type}{Type of transcript.}
}

\author{Gordon Smyth}

\seealso{
\code{\link{catchSalmon}}, \code{\link[limma]{strsplit2}}.
}

\concept{Transcript expression}
\concept{Reading data files}
