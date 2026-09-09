\name{catchSalmon}
\alias{catchSalmon}
\alias{catchKallisto}
\alias{catchOarfish}
\alias{catchRSEM}

\title{Process kallisto, Salmon, Oarfish, or RSEM Transcript Output}

\description{Read transcript counts from RSEM, kallisto, Salmon or Oarfish output for a series of biological samples and use bootstrap or posterior distribution samples to estimate the read-to-transcript-ambiguity for each transcript.}

\usage{
catchKallisto(parent.dir = NULL, sample.dirs = NULL,
          DGEList = TRUE, divide = FALSE,
          impute.eff.len = TRUE, offset.prior = TRUE, verbose = TRUE)
catchSalmon(parent.dir = NULL, sample.dirs = NULL,
          DGEList = TRUE, divide = FALSE,
          impute.eff.len = TRUE, offset.prior = TRUE, verbose = TRUE)
catchOarfish(parent.dir = NULL, prefixes = NULL,
          DGEList = TRUE, divide = FALSE, verbose = TRUE)
catchRSEM(parent.dir = NULL, files = NULL, ngibbs = 100,
          DGEList = TRUE, divide = FALSE, offset.prior = TRUE, verbose = TRUE)
}

\arguments{
\item{parent.dir}{name of parent directory containing the output directories from kallisto or Salmon or the output files from Oarfish or RSEM.}
\item{sample.dirs}{
  character vector giving names of or paths to the sample directories created by kallisto or Salmon.
  Each entry corresponds to one RNA-seq sample.
  If \code{parent.dir} and \code{sample.dirs} are both provided, then the \code{sample.dirs} are relative to \code{parent.dir}.
  If \code{sample.dirs} is \code{NULL}, then all the kallisto or Salmon directories found in \code{parent.dir} will be read.
}
\item{prefixes}{
  character vector giving the sample-specific file prefixes created by Oarfish.
  Each entry corresponds to one RNA-seq sample.
  Can contain a directory separator.
  By default, all files in the \code{parent.dir} directory with names ending in \code{.quant}, \code{.infreps.pq} and \code{.meta_info.json} will be read.
}
\item{files}{
  character vector specifying the \code{isoforms.results} files output by RSEM.
  Each entry corresponds to one RNA-seq sample.
  By default, all files in the \code{parent.dir} directory with names ending in \code{isoforms.results} will be read.
}
\item{DGEList}{logical. If \code{TRUE}, the function will output a \code{DGEList} instead of a plain list.}
\item{divide}{logical. If \code{TRUE}, counts are divided by the RTA overdispersions.}
\item{impute.eff.len}{
  whether to impute missing effective transcript lengths from Salmon.
}
\item{offset.prior}{logical. If \code{TRUE}, a prior offset matrix representing differences in effective transcript length will be included in the DGEList.}
\item{verbose}{logical. If \code{TRUE}, progress information will be sent to standard output as each sample is processed.}
\item{ngibbs}{number of Gibbs samples used to generate the RSEM results files.}
}

\details{
\emph{Reading transcript counts}.
\code{catchSalmon} assumes that Salmon (Patro et al 2017; Zakeri et al 2017) has been run to estimate transcript counts for one or more RNA samples and that either bootstrap resamples or Gibbs posterior samples are included in the Salmon output.
\code{catchSalmon} will automatically detect whether bootstrap or Gibbs samples are available.
The number of technical resamples does not need to be the same for each RNA sample, but the resamples should always be of the same type, either bootstrap or Gibbs.
We recommend that at least 200 bootstrap or Gibbs resamples are generated in total across all the RNA samples (Baldoni et al 2024b).

\code{catchKallisto} assumes that kallisto (Bray et al 2016) has been run to estimate transcript counts for one or more RNA samples and that bootstrap samples have also been generated.
The number of bootstrap resamples does not need to be the same for each RNA sample.
\code{catchKallisto} uses the rhdf5 package to read the HDF5 format files output by kallisto.

\code{catchRSEM} reads \code{*.isoforms.results} files created by RSEM (Li and Dewey, 2011), which should contain means and standard deviations from Gibbs posterior samples as well as the estimated transcript counts.
\code{catchRSEM} cannot detect the number of Gibbs samples performed for each RNA sample, so this must be specified via the \code{ngibbs} argument.
The number of Gibbs samples is assumed to be the same for each RNA sample.

\code{catchOarfish} assumes that Oarfish (Jousheghani et al 2025) has been run to estimate transcript counts for one or more RNA samples and that bootstrap resamples are included in the output.

\emph{RTA overdispersion}.
These functions read the transcript counts and use the technical (bootstrap or Gibbs) resamples to estimate an overdispersion parameter for each transcript.
The overdispersion represents the variance inflation that occurs from ambiguity in assigning sequence reads to transcripts, a phenomenon that is called \emph{read to transcript ambiguity} (RTA) overdispersion by Baldoni et al (2024ab).
Transcripts that overlap other transcripts and have greater quantification uncertainty will have larger RTA overdispersions.
The RTA overdispersions are greater than or equal to 1, with 1 representing no variance inflation.

\emph{Divided counts}.
To assess differential transcript expression, the transcript counts can be divided by the overdispersion parameters, after which the divided-counts can be input into standard differential expression pipelines designed for gene-level counts (Baldoni et al 2024a).
The edgeR quasi pipeline has been found to perform well with the divided counts (Baldoni et al 2024ab).
The divided-counts behave much like negative binomial counts and show the same mean-variance trends as for gene-level RNA-seq counts.
While bootstrap and Gibbs resamples both perform well, Baldoni et al (2024b) show that Gibbs resamples are computationally faster than bootstrap resamples and give slightly better performance in the downstream differential expression analyses.

\emph{Imputing missing transcript lengths}.
For paired-end RNA-seq, Salmon uses sample-specific empirical fragment distributions to estimate the effective transcript lengths.
For very short transcripts, there may be too few entries in the empirical distribution to estimate the conditional fragment length.
In such cases, Salmon falls back to setting the effective transcript length equal to the actual length even though the effective transcript length should in principle be much shorter.
If \code{impute.eff.len} is \code{TRUE}, then any effective transcript length equal to the transcript's actual length will be reset to the minimum effective transcript length for that transript across all the samples in the series.

\emph{Gencode imbedded annotation}.
Gencode packs transcript and gene annotation into composite transcript IDs, with fields separated by vertical bar \code{|} characters.
All the functions above will detect Gencode-style transcript names and will unpack the annotation fields into separate annotation columns using \code{splitGencodeTxNames}.
}

\value{
If \code{DGEList=TRUE}, a \code{\link[edgeR:DGEList-class]{DGEList}} object with components:
\item{counts}{matrix of transcript counts, with rows for transcripts and columns for RNA samples.}
\item{genes}{data.frame of transcript information with columns \code{Length}, \code{AveEffLen}, \code{Max2MinEffLen}, and \code{Overdispersion} corresponding to transcript length, effective transcript length averaged over samples, ratio of largest to smallest length for different samples, and RTA overdispersion respectively.}
\item{overdispersion.prior}{median overdispersion, used to moderate the transcript-wise overdispersion values.}
\item{resample.type}{character vector giving type of resampling (\code{"bootstrap"} or \code{"gibbs"}) for each sample.}
\item{divided.counts}{a logical indicating whether the counts have been divided by the RTA overdispersions.}
\item{offset.prior}{matrix of prior offsets representing relative differences in effective lengths for each transcript on the log-scale (only if \code{offset.prior=TRUE}). Rows add to zero.}
\item{other$effective.length}{matrix of effective transcript lengths.}

If Gencode-style composite rownames are detected, then the rownames will be reset to plain Ensembl transcript IDs, without veresion number suffixes, and the \code{genes} data.frame will include columns:
  \item{EnsG}{Ensembl gene ID, without version number suffix.}
  \item{HavG}{Havana gene ID.}
  \item{HavT}{Havana transcript ID.}
  \item{TxName}{Transcript name.}
  \item{GeneName}{Gene name.}
  \item{Type}{Type of transcript.}

If \code{DGEList=FALSE}, then a plain list is created.
The components are as documented above for DGELists except that the annotation data.frame is called \code{annotation} instead of \code{genes},
and the \code{effective.length} matrix is a top-level component.

For \code{catchOarfish}, the \code{annotation} data.frame does not contain the \code{AveEffLen} and \code{Max2MinEffLen} columns.
}

\author{Gordon Smyth and Pedro Baldoni}

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
\doi{10.1093/nar/gkaf1305}.

Bray NL, Pimentel H, Melsted P, Pachter L (2016).
Near-optimal probabilistic RNA-seq quantification.
\emph{Nature Biotechnology} 34(5), 525-527.

Jousheghani ZZ, Singh NP, Patro R (2025).
Oarfish: enhanced probabilistic modeling leads to improved accuracy in long read transcriptome quantification.
\emph{Bioinformatics} 41, i304-i313.
\doi{10.1093/bioinformatics/btaf240}.

Li B, Dewey CN (2011).
RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome.
\emph{BMC Bioinformatics} 12, 323.

Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C (2017).
Salmon provides fast and bias-aware quantification of transcript expression.
\emph{Nature Methods} 14(4), 417-419.

Zakeri M, Srivastava A, Almodaresi F, Patro R (2017).
Improved data-driven likelihood factorizations for transcript abundance estimation.
\emph{Bioinformatics} 33(14), i142-i151.
}

\seealso{
See \code{\link{catchSalmonGene}} for gene-level counts and overdispersions.
}

\examples{
\dontrun{
# Read Salmon ouput and estimate overdispersion for each transcript
s <- catchSalmon(paths)

# Divide the transcript counts by the overdispersions ready for a standard edgeR DE analysis
dge <- DGEList(counts=s$counts/s$annotation$Overdispersion, genes=s$annotation)
}}

\concept{Transcript expression}
\concept{Reading data files}
