\name{readBismark2PC}
\alias{readBismark2PC}

\title{Read Bismark Coverage Files}

\description{Read Bismark coverage files containing methylated and unmethylated read counts for CpG loci and create a PCList.}

\usage{readBismark2PC(files, sample.names=NULL, readr=TRUE, verbose=TRUE)}

\arguments{ 
\item{files}{character vector of file names.}
\item{sample.names}{character vector of sample names. If \code{NULL}, sample names will be extracted from the file names.}
\item{readr}{logical. If \code{TRUE}, \code{readr} package is used to read the coverage files, otherwise \code{read.delim} is used.}
\item{verbose}{logical. If \code{TRUE}, read progress messages are sent to standard output.}
}

\details{
This function reads tab-delimited coverage files output by Bismark software.
Counts from multiple files are collated into a \code{\link[edgeR:PCList-class]{PCList}} object.
}

\note{
This function represents genomic loci as integers, so the largest locus position must be less than the maximum integer in R (about \code{2e9}).
The number of chromosomes times the largest locus position must be less than \code{1e16}.
}

\value{
A \code{\link[edgeR:PCList-class]{PCList}} object with a row for each unique genomic locus found in the files and two count matrices (containing methylated and unmethylated counts) for each sample.
}

\author{Gordon Smyth}

\references{
Chen Y, Pal B, Visvader JE, Smyth GK (2017).
Differential methylation analysis of reduced representation bisulfite sequencing experiments using edgeR.
\emph{F1000Research} 6, 2055.
\doi{10.12688/f1000research.13196.2}.
}

\examples{
# Two small Bismark coverage files: chr, start, end, meth percentage, count_meth, count_unmeth
tab <- intToUtf8(9)
f1 <- tempfile(fileext=".cov")
f2 <- tempfile(fileext=".cov")
writeLines(paste("chr1", c(10,20), c(10,20), c(50,80), c(5,8), c(5,2), sep=tab), f1)
writeLines(paste("chr1", c(10,30), c(10,30), c(60,70), c(6,7), c(4,3), sep=tab), f2)
y <- readBismark2PC(c(f1,f2), readr=FALSE, verbose=FALSE)
dim(y)
y$counts
}

\concept{Reading data files}
\concept{Differential methylation}
