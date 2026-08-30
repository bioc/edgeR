\name{tpmProp}
\alias{tpmProp}

\title{TPM Proportion by Transcript for Each Gene}

\description{Compute the expression proportion from each transcript for each gene.}

\usage{
tpmProp(tpm, geneid)
}

\arguments{ 
\item{tpm}{numeric matrix of TPM values, with rows corresponding to transcripts.}
\item{geneid}{character vector or factor of Gene IDs.}
}

\details{
Convert a matrix of TPMs to expression proportions for each gene.
}

\value{
A numeric matrix of TPM proportions, adding to 1 for each gene in each sample.
}

\references{
Baldoni PL, Chen L, Li M, Chen Y, Smyth GK (2025).
Dividing out quantification uncertainty enables assessment of differential transcript usage with limma and edgeR.
\emph{bioRxiv}
\doi{10.1101/2025.04.07.647659}.
}

\author{Gordon Smyth}

\seealso{
\code{\link{tpm}}
}

\keyword{rna-seq}
\concept{Differential transcript usage}
\concept{Differential usage}
