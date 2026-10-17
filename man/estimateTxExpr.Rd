\name{estimateTxExpr}
\alias{estimateTxExpr}

\title{Estimate Transcript Expression Levels from Equivalence Class Counts}

\description{Estimate transcript levels given transcript lengths, equivalence class counts, and compatibility matrix.}

\usage{
estimateTxExpr(lengths, counts, compatibility, maxit = 1000L, eps = 1e-6, start.rpb = NULL)
}

\arguments{
\item{lengths}{vector of effective transcript lengths, one for each transcript.}
\item{counts}{
  vector of equivalence class read counts, one for each equivalence class.
  Alternatively, can be a matrix where each column is a vector of equivalence class reads.
  In this case, the function will return transcript expression estimates for each input column.
}
\item{compatibility}{list of same length as \code{counts}. Each entry is a vector of indices of transcripts belonging to the equivalence class, i.e., transcript compatibility with these reads.}
\item{maxit}{maximum number of EM algorithm iterations allowed.}
\item{eps}{stopping criterion.}
\item{start.rpb}{starting values for the tx expression values.}
}

\details{
Estimate transcript expression levels by the EM algorithm from equivalence class counts.
Expression levels are in reads per base, proportional to TPM for each transcript.
In principle, the quantifications should be similar to those returned by kallisto for an inter-connected set of transcripts.

The \code{counts} argument can be a matrix, to facilitate simulations in which each column is a simulated vector of equivalence class counts.
}

\value{
A list the following components
\item{rpb}{vector of transcript expression levels in units of reads per base, or a matrix with \code{ncol(counts)} columns.}
\item{est.count}{vector of estimated transcript counts, or a matrix with \code{ncol(counts)} columns.}
\item{niter}{number of EM iterations used.}
\item{delta}{relative accuracy achieved.}
\item{contraction.factor}{rate of convergence of the EM algorithm. The maximum difference between the working estimate and the maximum likelihood estimate contracts by this factor each iteration (Smyth 1996).}
}

\author{Gordon Smyth}

\references{
Dempster AP, Laird NM, Rubin DB (1977).
Maximum likelihood from incomplete data via the EM algorithm.
\emph{Journal of the Royal Statistical Society Series B} 39(1), 1-22.

Smyth GK (1996).
Partitioned algorithms for maximum likelihood and other nonlinear estimation.
\emph{Statistics and Computing} 6(3), 201-216.
\doi{10.1007/BF00140865}.
}

\examples{
# Two overlapping transcripts, one long and one short.
lengths <- c(Tx1=2000,Tx2=500)
counts <- c(EC1=100,EC2=50)
compatibility <- list(EC1=1,EC2=c(1,2))
estimateTxExpr(lengths, counts, compatibility)

# Multiple input count vectors
counts <- cbind(c(100,50),c(110,40))
rownames(counts) <- c("EC1","EC2")
colnames(counts) <- c("S1","S2")
estimateTxExpr(lengths, counts, compatibility)
}

\concept{Transcript expression}
