\name{catchSalmon}
\alias{catchSalmon}
\alias{catchKallisto}
\alias{catchRSEM}

\title{Process Salmon, kallisto or RSEM Transcript Output}

\description{Read transcript counts from RSEM, kallisto or Salmon output for a series of biological samples and use bootstrap or posterior distribution samples to estimate the read-to-transcript-ambiguity for each transcript.}

\usage{
catchSalmon(paths, verbose = TRUE)
catchKallisto(paths, verbose = TRUE)
catchRSEM(files = NULL, ngibbs = 100, verbose = TRUE)
}

\arguments{
\item{paths}{character vector giving paths to the directories created by a pseudoaligner (either kallisto or Salmon). Each entry corresponds to one RNA-seq sample.}
\item{verbose}{logical. If \code{TRUE}, progress information will be sent to standard output as each sample is processed.}
\item{files}{character vector specifying the \code{isoforms.results} files output by RSEM. Each entry corresponds to one RNA-seq sample. By default, all files in the current working directory with names ending in \code{isoforms.results} will be read.}
\item{ngibbs}{number of Gibbs samples used to generate the RSEM results files.}
}

\details{
\code{catchSalmon} assumes that Salmon (Patro et al 2017; Zakeri et al 2017) has been run to estimate transcript counts for one or more RNA samples and that either bootstrap resamples or Gibbs posterior samples are included in the Salmon output.
\code{catchSalmon} will automatically detect whether bootstrap or Gibbs samples are available.
The number of technical resamples does not need to be the same for each RNA sample, but the resamples should always be of the same type, either bootstrap or Gibbs.
We recommend that at least 200 bootstrap or Gibbs resamples are generated in total across all the RNA samples (Baldoni et al 2024b).

\code{catchKallisto} assumes that kallisto (Bray et al 2016) has been run to estimate transcript counts for one or more RNA samples and that bootstrap samples have also been generated.
The number of bootstrap resamples does not need to be the same for each RNA sample.

\code{catchRSEM} reads \code{*.isoforms.results} files created by RSEM (Li and Dewey, 2011), which should contain means and standard deviations from Gibbs posterior samples as well as the estimated transcript counts.
\code{catchRSEM} cannot detect the number of Gibbs samples performed for each RNA sample, so this must be specified via the \code{ngibbs} argument.
The number of Gibbs samples is assumed to be the same for each RNA sample.

These functions read the transcript counts and use the technical (bootstrap or Gibbs) resamples to estimate an overdispersion parameter for each transcript.
The overdispersion represents the variance inflation that occurs from ambiguity in assigning sequence reads to transcripts, a phenomenon that is called \emph{read to transcript ambiguity} (RTA) overdispersion by Baldoni et al (2024ab).
Transcripts that overlap other transcripts and have greater quantification uncertainty will have larger RTA overdispersions.
The RTA overdispersions are greater than or equal to 1, with 1 representing no variance inflation.

To assess differential transcript expression, the transcript counts can be divided by the overdisperson parameters, after which the scaled counts can be input into standard differential expression pipelines designed for gene-level counts (Baldoni et al 2024a).
The edgeR quasi pipeline has been found to perform well with the scaled counts (Baldoni et al 2024ab).
The scaled counts behave much like negative binomial counts and show the same mean-variance trends as for gene-level RNA-seq counts.
While bootstrap and Gibbs resamples both perform well, Baldoni et al (2024b) show that Gibbs resamples are computationally faster than bootstrap resamples and give slightly better performance in the downstream differential expression analyses.
}

\value{
A list containing components
\item{counts}{matrix of transcript counts, with rows for transcripts and columns for RNA samples.}
\item{annotation}{data.frame of transcript information with columns \code{Length}, \code{EffectiveLength}, and \code{Overdispersion} corresponding to transcript length, effective transcript length and RTA overdispersion respectively.}
\item{overdispersion.prior}{median overdispersion, used to moderate the transcript-wise overdispersion values.}
\item{resample.type}{character vector giving type of resampling (\code{"bootstrap"} or \code{"gibbs"}) for each sample. Only for \code{catchSalmon}.}
}

\author{Gordon Smyth and Pedro Baldoni}

\references{
Baldoni PL, Chen Y, Hediyeh-zadeh S, Liao Y, Dong X, Ritchie ME, Shi W, Smyth GK (2024a).
Dividing out quantification uncertainty allows efficient assessment of differential transcript expression with edgeR.
\emph{Nucleic Acids Research} 52(3), e13.
\doi{10.1093/nar/gkad1167}.

Baldoni PL, Chen L, Smyth GK (2024b).
Faster and more accurate assessment of differential transcript expression with Gibbs sampling and edgeR 4.0.
\emph{bioRxiv}
\doi{10.1101/2024.06.25.600555}.

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

