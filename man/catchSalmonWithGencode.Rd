\name{catchSalmonWithGencode}
\alias{catchSalmonWithGencode}

\title{Process Salmon Output With Gencode Annotation to Produce Gene Counts}

\description{Produce a gene count object from Salmon output for a series of biological samples with Gencode annotation.}

\usage{
catchSalmonWithGencode(parent.dir = NULL, sample.dirs = NULL,
        DGEList = TRUE, divide = FALSE, impute.eff.len = TRUE,
        offset.prior = TRUE, gene.length = "moderate", verbose = TRUE)
}

\arguments{
\item{parent.dir}{name of parent directory containing the Salmon output directories.}
\item{sample.dirs}{
  character vector giving paths to the sample directories created by Salmon.
  If \code{parent.dir} and \code{sample.dirs} are both provided, then the \code{sample.dirs} are relative to \code{parent.dir}.
  If \code{sample.dirs} is \code{NULL}, then all the Salmon directories found in \code{parent.dir} will be read.
}
\item{DGEList}{logical. If \code{TRUE}, the function will output a \code{DGEList} instead of a plain list.}
\item{divide}{logical. If \code{TRUE}, counts are divided by the RTA overdispersions.}
\item{impute.eff.len}{
  whether to impute missing effective transcript lengths from Salmon.
}
\item{offset.prior}{logical. If \code{TRUE}, a prior offset matrix representing differences in effective transcript length will be included in the DGEList.}
\item{gene.length}{How to compute average gene length. Possibilities are \code{"moderate"}, \code{"tximport"} and \code{"simple"}.}
\item{verbose}{logical. If \code{TRUE}, progress information will be sent to standard output as each sample is processed.}
}

\details{
This function assumes that Salmon (Patro et al 2017; Zakeri et al 2017) has been run with a Gencode transcript FASTA file to estimate transcript counts for one or more RNA samples, and that either bootstrap resamples or Gibbs posterior samples are included in the Salmon output.
The function will automatically detect whether bootstrap or Gibbs samples are available.

The function will unpack the Gencode transcript names using \code{splitGencodeTxNames} to obtain gene and transcript IDs.
The transcript-level counts will be collapsed to gene-level counts and RTA overdispersions will be estimated at the gene level.

Average transcript lengths are computed for each gene and optionally used to create a prior offset matrix.
The average lengths are computed by three possible methods as specified by \code{gene.length}, either a simple average of the transcript lengths (\code{"simple"}), or the method of Soneson et al (2016) (\code{"tximport"}),
or by a method that moderates the lengths weakly towards the gene average and, even more weakly, towards the simple average of transcript lengths (\code{"moderate"}).
The moderated method is recommended.

For paired-end RNA-seq, Salmon uses sample-specific empirical fragment distributions to estimate the effective transcript lengths.
For very short transcripts, there may be too few entries in the empirical distribution to estimate the conditional fragment length.
In such cases, Salmon falls back to setting the effective transcript length equal to the actual length even though the effective transcript length should in principle be much shorter.
If \code{impute.eff.len} is \code{TRUE}, then any effective transcript length equal to the transcript's actual length will be reset to the minimum effective transcript length for that transript across all the samples in the series.
}

\value{
If \code{DGEList=FALSE}, a plain list containing components
\item{counts}{matrix of transcript counts, with rows for transcripts and columns for RNA samples.}
\item{length}{matrix of effective transcript lengths, with rows for transcripts and columns for RNA samples.}
\item{annotation}{data.frame of transcript information with columns \code{Length}, \code{AveLength}, \code{Max2MinLength}, and \code{Overdispersion} corresponding to transcript length, effective transcript length averaged over samples, ratio of largest to smallest length for different samples, and RTA overdispersion respectively.}
\item{overdispersion.prior}{median overdispersion, used to moderate the transcript-wise overdispersion values.}
\item{resample.type}{character vector giving type of resampling (\code{"bootstrap"} or \code{"gibbs"}) for each sample.}
\item{divided.counts}{a logical indicating whether the counts have been divided by the RTA overdispersions.}

If \code{DGEList=TRUE}, a \code{\link[edgeR:DGEList-class]{DGEList}} object as the output.
In this case, the \code{genes} component of the object will contain the same columns as in the \code{annotation} data.frame from the plain list output.
If \code{offset.prior=TRUE}, then the \code{offset.prior} component will contain prior offsets representing sample-wise differences in effective transcript lengths.
As well as standard DGEList components, the output DGEList will contain components \code{overdispersion.prior}, \code{resample.type}, and \code{divided.counts} equivalent to those from the plain list output.

The rownames will be reset to plain Ensembl gene IDs and the \code{genes} data.frame will include columns:
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

Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C (2017).
Salmon provides fast and bias-aware quantification of transcript expression.
\emph{Nature Methods} 14(4), 417-419.

Soneson C, Love MI, Robinson MD (2016).
Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences.
\emph{F1000Research} 4, 1521.
\doi{10.12688/f1000research.7563.2}.

Zakeri M, Srivastava A, Almodaresi F, Patro R (2017).
Improved data-driven likelihood factorizations for transcript abundance estimation.
\emph{Bioinformatics} 33(14), i142-i151.
}

\seealso{
\code{\link{catchSalmon}}, \code{\link{splitGencodeTxNames}}.
}

\examples{
\dontrun{dge <- DGEList(mydir)}
}

\concept{Transcript expression}
\concept{Reading data files}
