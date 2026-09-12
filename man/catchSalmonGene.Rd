\name{catchSalmonGene}
\alias{catchKallistoGene}
\alias{catchSalmonGene}

\title{Process Kallisto or Salmon Output to Produce Gene Counts}

\description{Produce a gene count object from kallisto or Salmon output for a series of biological samples, optionally with Gencode annotation.}

\usage{
catchKallistoGene(parent.dir = NULL, sample.dirs = NULL,
        tx2gene = NULL, remove.version.numbers = TRUE,
        DGEList = TRUE, divide = FALSE, impute.eff.len = TRUE,
        offset.prior = TRUE, gene.length = "moderate", verbose = TRUE)
catchSalmonGene(parent.dir = NULL, sample.dirs = NULL,
        tx2gene = NULL, remove.version.numbers = TRUE,
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
\item{tx2gene}{data.frame mapping transcript IDs to gene IDs, with first column for transcript IDs and the second for the corresponding gene IDs. Not needed if Gencode annotation was used for the quantifications.}
\item{remove.version.numbers}{logical, whether to remove version numbers from transcript and gene IDs.}
\item{DGEList}{logical, whether to output a \code{DGEList} instead of a plain list.}
\item{divide}{logical, whether to divide the gene counts by the RTA overdispersions.}
\item{impute.eff.len}{logical, whether to impute missing effective transcript lengths.}
\item{offset.prior}{logical, whether to include a prior offset matrix representing differences in effective gene length.}
\item{gene.length}{how to compute effective gene length. Possibilities are \code{"moderate"}, \code{"tximport"} and \code{"simple"}.}
\item{verbose}{logical, whether to message progress information to standard output as each sample is processed.}
}

\details{
\code{catchSalmonGene()} assumes that Salmon (Patro et al 2017; Zakeri et al 2017) has been run to estimate transcript counts for one or more RNA samples, and that either bootstrap resamples or Gibbs posterior samples are included in the Salmon output.
The function reads the Salmon quantifications and summarizes them at the gene level.
The function will automatically detect whether bootstrap or Gibbs samples are available.

\code{catchKallistoGene()} assumes that kallisto (Bray et al 2016) has been run to estimate transcript counts for one or more RNA samples and that bootstrap samples have also been generated.
The number of bootstrap resamples does not need to be the same for each RNA sample.
\code{catchKallistoGene()} uses the rhdf5 package to read the HDF5 format files output by kallisto.

The transcript-level counts will be collapsed to gene-level counts and RTA overdispersions will be estimated at the gene level, using a method similar at Baldoni et al (2024ab) but applied to gene level counts instead of transcript counts.

The function needs to know which gene each transcript belongs to.
If Salmon was run with a Gencode transcript FASTA file, then \code{catchSalmonGene} will detect that automatically.
It will then unpack the Gencode transcript names using \code{splitGencodeTxNames} to obtain gene and transcript IDs.
Otherwise, a data.frame \code{tx2gene} needs to be provided to map transcript to gene IDs.

Average transcript lengths are computed for each gene and optionally used to create a prior offset matrix.
The average lengths are computed by three possible methods as specified by \code{gene.length}, either a simple average of the transcript lengths (\code{"simple"}), or the method of Soneson et al (2016) (\code{"tximport"}),
or by a method that moderates the lengths weakly towards the gene average and, even more weakly, towards the simple average of transcript lengths (\code{"moderate"}).
The moderated method is recommended.

For paired-end RNA-seq, kallisto and Salmon use sample-specific empirical fragment distributions to estimate the effective transcript lengths.
For very short transcripts, there may be too few entries in the empirical distribution to estimate the conditional fragment length.
In such cases, kallisto and Salmon fall back to setting the effective transcript length equal to the actual length even though the effective transcript length should in principle be much shorter.
If \code{impute.eff.len} is \code{TRUE}, then any effective transcript length equal to the transcript's actual length will be reset to the minimum effective transcript length for that transcript across all the samples in the series.
}

\value{
If \code{DGEList=TRUE}, a \code{\link[edgeR:DGEList-class]{DGEList}} object with gene IDs as row.names and the following components:
\item{counts}{matrix of gene counts, with rows for genes and columns for RNA samples.}
\item{genes}{data.frame of gene annotation.}
\item{genes$GeneName}{gene name (if Gencode annotation detected).}
\item{genes$NTx}{number of transcripts per gene.}
\item{genes$MaxTxLen}{maximum transcript length per gene.}
\item{genes$AveEffLen}{effective gene length, averaged over samples.}
\item{genes$Max2MinEffLen}{ratio of largest to smallest effective gene length across samples.}
\item{genes$Overdispersion}{genewise RTA overdispersion.}
\item{overdispersion.prior}{median of genewise overdispersions, used to moderate the overdispersion values.}
\item{resample.type}{character vector giving the resampling type (\code{"bootstrap"} or \code{"gibbs"}) for each sample.}
\item{divided.counts}{logical, indicating whether the counts have been divided by the RTA overdispersions.}
\item{offset.prior}{matrix of prior offsets representing relative differences in effective lengths for each gene on the log-scale (only if \code{offset.prior=TRUE}). Rows add to zero.}
\item{other$effective.length}{matrix of effective gene lengths.}
\item{other$tpm}{matrix of genewise TPM values.}

If \code{DGEList=FALSE}, then a plain list is created.
The components are as documented above for DGELists except that the annotation data.frame is called \code{annotation} instead of \code{genes}, and the \code{effective.length} and \code{tpm} matrices are top-level components.
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
\code{\link{catchKallisto}}, \code{\link{catchSalmon}}, \code{\link{splitGencodeTxNames}}.
}

\examples{
\dontrun{dge <- DGEList(mydir)}
}

\concept{Transcript expression}
\concept{Reading data files}
