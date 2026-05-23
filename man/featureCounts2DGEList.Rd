\name{featureCounts2DGEList}
\alias{featureCounts2DGEList}

\title{Convert featureCounts object to a DGEList}

\description{
Converts the list output by \code{Rsubread::featureCounts} to a DGEList object.
}

\usage{
featureCounts2DGEList(x, verbose = TRUE)
}

\arguments{
  \item{x}{a list produced by \code{Rsubread::featureCounts.}}
  \item{verbose}{if \code{TRUE}, will give a message to report on the number of junctions that can be assigned to a gene.}
}

\details{
Rsubread's \code{featureCounts} function counts reads by features (typically exons or genomic intervals) or meta-features (typically genes).
It can also count reads crossing exon-exon junctions or reads internal to exons.
This function assembles all the read counts output by \code{featureCounts} into an DGEList, ensuring unique row.names where appropriate.
The proportion of reads assigned to features is also stored.

If \code{x} contains exon and junction counts, then exons and junctions are numbered within genes, and rownames are of the form \code{GeneID.e1} for exons and \code{GeneID.j1} for junctions.
Junctions are assigned length 1 in the output.

This function is designed for gene-orientated analyses, including differential usage analyses of exons and junctions.
If \code{x} contains junction reads, then only junctions for which both the donor and acceptor sites can be assigned to the same gene are included in the output object.
The function also drops the columns \code{Transcript} and \code{Status} from \code{x$counts_junction}.
Users who want to retain all junctions and annotation, including chimeric junctions, can use \code{DGEList(x$counts_junction)}.
}

\value{
A \code{\link[edgeR:DGEList-class]{DGEList}} data object.
The \code{samples} component of the object includes a column \code{PropAssigned} recording the proportion of reads (or read-pairs) successfully assigned to features.
The \code{genes} component of the object includes columns \code{GeneID}, \code{Chr}, \code{Start}, \code{End}, \code{Strand}, \code{Length}, and \code{Exon}.
The \code{Exon} column gives the numbering of the exons and junctions within genes.
}

\author{Gordon Smyth.}

\seealso{
\code{\link[edgeR:DGEList-class]{DGEList-class}} defines DGEList objects.

See the \code{featureCounts} function in the Rsubread package for creation of \code{x}.
}

\concept{edgeR classes}
\concept{differential transcript usage}
\concept{Differential usage}
