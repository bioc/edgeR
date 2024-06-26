\name{catchSalmon}
\alias{catchSalmon}
\alias{catchKallisto}
\alias{catchRSEM}

\title{Process Salmon, kallisto or RSEM Transcript Output}

\description{Read transcript counts from kallisto or Salmon output for a series of biological samples and use bootstrap samples to estimate the read-to-transcript ambiguity for each transcript.}

\usage{
catchKallisto(paths, verbose = TRUE)
catchSalmon(paths, verbose = TRUE)
catchRSEM(files = NULL, ngibbs = 100, verbose = TRUE)
}

\arguments{
\item{paths}{character vector giving paths to the directories created by the pseudoaligner (either kallisto or Salmon). Each entry corresponds to one RNA-seq sample.}
\item{verbose}{logical. If \code{TRUE}, summary information is shown as each sample is processed.}
\item{files}{character vector specifying the isoforms.results files output by RSEM. Each entry corresponds to one RNA-seq sample. By default, all files in the current directory with names ending in \code{isoforms.results} will be read.}
\item{ngibbs}{number of Gibbs samples used to generate the RSEM results files.}
}

\details{
The functions \code{catchSalmon} or \code{catchKallisto} assume that kallisto (Bray et al 2016) or Salmon (Patro et al 2017; Zakeri et al 2017) have been run to obtain estimated transcript counts for one or more RNA samples and that bootstrap samples or Gibbs samples have also been generated.
These functions catch the counts and use the bootstrap samples to estimate an overdispersion parameter for each transcript.
The overdispersion represents the variance inflation that occurs from ambiguity in assigning sequence reads to transcripts.
Transcripts that overlap other transcripts and have greater quantification uncertainty will have larger overdispersions.

\code{catchRSEM} reads \code{*.isoforms.results} files created by RSEM (Li and Dewey, 2011), which should contain means and standard deviations from the Gibbs posterior samples as well as the estimated transcript counts.

To assess differential transcript expression, the transcript counts can be divided by the overdisperson parameters, after which the scaled counts can be input into standard differential expression pipelines designed for gene-level counts (Baldoni et al 2023).
The edgeR quasi pipeline has been found to perform well with the scaled counts (Baldoni et al 2024).
The scaled counts behave much like negative binomial counts and show the same mean-variance trends as for gene-level RNA-seq counts.
}

\value{
A list containing components
\item{counts}{matrix of transcript counts}
\item{annotation}{data.frame of transcript information with columns \code{Length}, \code{EffectiveLength} and \code{Overdispersion}.}
\item{overdispersion.prior}{median overdispersion, used to moderate the transcript-wise overdispersion values.}
\item{resample.type}{character vector giving type of resampling (\code{"bootstrap"} or \code{"gibbs"}) for each sample. Only for \code{catchSalmon}.}
}

\author{Gordon Smyth and Pedro Baldoni}

\references{
Baldoni PL, Chen Y, Hediyeh-zadeh S, Liao Y, Dong X, Ritchie ME, Shi W, Smyth GK (2024).
Dividing out quantification uncertainty allows efficient assessment of differential transcript expression with edgeR.
\emph{Nucleic Acids Research} 52(3), e13.
\doi{10.1093/nar/gkad1167}.

Bray NL, Pimentel H, Melsted P, Pachter L (2016).
Near-optimal probabilistic RNA-seq quantification.
\emph{Nature Biotechnology}, 34(5), 525-527.

Li B, Dewey CN (2011).
RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome.
\emph{BMC bioinformatics}, 12, 323.

Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C (2017).
Salmon provides fast and bias-aware quantification of transcript expression.
\emph{Nature Methods}, 14(4), 417-419.

Zakeri M, Srivastava A, Almodaresi F, Patro R (2017).
Improved data-driven likelihood factorizations for transcript abundance estimation.
\emph{Bioinformatics} 33(14), i142-i151.
}

\examples{
\dontrun{
# Read Salmon ouput and estimate overdispersion for each transcript
s <- catchSalmon(paths)

# Scale the transcript counts ready for a standard edgeR DE analysis
dge <- DGEList(counts=s$counts/s$annotation$Overdispersion, genes=s$annotation)
}}

\concept{Transcript expression}

