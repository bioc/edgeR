\name{catchSalmon}
\alias{catchSalmon}
\alias{catchKallisto}

\title{Process Kallisto or Salmon Output}

\description{Read transcriptwise counts from kallisto or Salmon output for a series of samples and use the bootstrap samples to estimate the mapping uncertainty for each transcript.}

\usage{
catchKallisto(paths, verbose = TRUE)
catchSalmon(paths, verbose = TRUE)
}

\arguments{
\item{paths}{character vector giving paths to the directories created by kallisto.}
\item{verbose}{logical. If \code{TRUE}, summary information is printed as each sample is catch.}
}

\details{
These functions assume that kallisto or Salmon have been run to obtain estimated transcript counts for one or more RNA samples, and that bootstrap samples have also been generated.
These functions catch the counts and use the bootstrap samples to estimate an over-dispersion coefficient for each transcript.
Transcripts that overlap other transcripts and have greater read mapping uncertaintly will have greater over-dispersion coefficients.

The data is then ready for analysis in edgeR.
}

\value{
A list containing components
\item{counts}{matrix of transcript counts}
\item{annotation}{data.frame of transcript information.}
}

\author{Gordon Smyth}

\references{
Baldoni PL, Chen Y, Hediyeh-zadeh S, Liao Y, Dong X, Ritchie ME, Shi W, Smyth GK (2023).
Dividing out quantification uncertainty allows efficient assessment of differential transcript expression.
\emph{bioRxiv}
\doi{10.1101/2023.04.02.535231}.
}

\examples{
\dontrun{
s <- catchSalmon(paths)
dge <- DGEList(counts=s$counts/s$annotation$Overdispersion, genes=s$annotation)
}}

\concept{Transcript expression}

