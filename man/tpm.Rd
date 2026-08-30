\name{tpm}
\alias{tpm}

\title{Transcripts per Million}

\description{Compute transcripts per million (TPM) from counts or from a fitted model object.}

\usage{
tpm(y, effective.tx.length, rta.overdispersion=NULL, shrunk=FALSE, ...)
}

\arguments{ 
\item{y}{
  any object that can be input to \code{cpm} with rows corresponding to transcripts.
  Typically a DGEList containing divided transcript counts or a DGEGLM fitted model object fitted to divided transcript counts.
}
\item{effective.tx.length}{numeric vector giving effective transcript lengths.}
\item{rta.overdispersion}{numeric vector giving overdispersions used to compute divided counts.}
\item{shrunk}{logical, if \code{y} is a DGEGLM fitted model object, should shrunk or unshrunk coefficients be used? Default is unshrunk.}
\item{\dots}{other arguments are passed to \code{cpm}.}
}

\details{
This function computes transcripts per million (TPM) by rescaling counts-per-million at the transcript level by the effective transcript length and by any count scaling factor.
The function will typically be run as part of the DTU analysis after running \code{catchSalmon}, dividing the transcript counts by the RTA-overdispersion, and running \code{glmQLFit} on the divided counts.
The \code{tpm} function can then extract TPM values from the estimated model coefficients.

Note that the TPM values do not sum exactly to one million for each sample.
The total TPM varies from sample to sample because the effective library sizes estimated by \code{normLibSizes} imply that the total RNA production may differ from sample to sample.
The TPM values preserve the same fold-changes between samples and conditions as obtained from \code{cpm}.
}

\value{
A numeric matrix of TPM values.
}

\references{
Baldoni PL, Chen L, Li M, Chen Y, Smyth GK (2025).
Dividing out quantification uncertainty enables assessment of differential transcript usage with limma and edgeR.
\emph{bioRxiv}
\doi{10.1101/2025.04.07.647659}.
}

\author{Gordon Smyth}

\seealso{
\code{\link{cpm}}, \code{\link{diffSplice.DGEGLM}}
}

\keyword{rna-seq}
\concept{Differential transcript usage}
\concept{Differential usage}
