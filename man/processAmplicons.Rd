\name{processAmplicons}
\alias{processAmplicons}

\title{Process FASTQ files from pooled genetic sequencing screens}

\description{
Given a list of sample-specific index (barcode) sequences and hairpin/sgRNA-specific sequences from an amplicon sequencing screen, generate a DGEList of counts from the raw FASTQ files containing the sequence reads. The position of the index sequences and hairpin/sgRNA sequences is considered variable, with the hairpin/sgRNA sequences assumed to be located after the index sequences in the read.
}

\usage{
processAmplicons(readfile, readfile2 = NULL, barcodefile, hairpinfile,
                 allowMismatch = FALSE, barcodeMismatchBase = 1,
                 hairpinMismatchBase = 2, dualIndexForwardRead = FALSE,
                 verbose = FALSE, barcodesInHeader = FALSE,
                 hairpinBeforeBarcode = FALSE, plotPositions = FALSE)
}

\arguments{
\item{readfile}{character vector giving one or more FASTQ filenames. Must be uncompressed text files.}
\item{readfile2}{optional character vector giving one or more FASTQ filenames for reverse read. Must be uncompressed text files.}
\item{barcodefile}{filename containing sample-specific barcode IDs and sequences. File may be gzipped.}
\item{hairpinfile}{filename containing hairpin/sgRNA-specific IDs and sequences. File may be gzipped.}
\item{allowMismatch}{logical, indicates whether sequence mismatch is allowed.}
\item{barcodeMismatchBase}{maximum number of base sequence mismatches allowed in a barcode sequence when \code{allowMismatch} is \code{TRUE}.}
\item{hairpinMismatchBase}{maximum number of base sequence mismatches allowed in a hairpin/sgRNA sequence when \code{allowMismatch} is \code{TRUE}.}
\item{dualIndexForwardRead}{logical, indicates if forward reads contains a second barcode sequence (must be present in \code{barcodefile}) that should be matched.}
\item{verbose}{if \code{TRUE}, output program progress.}
\item{barcodesInHeader}{logical, indicates if barcode sequences should be matched in the header (sequence identifier) of each read (i.e. the first of every group of four lines in the FASTQ files).}
\item{hairpinBeforeBarcode}{logical, indicates that hairpin sequences appear before barcode sequences in each read. Setting this to TRUE will allow hairpins to exist anywhere in a read irrespective of barcode position. Leaving as FALSE will improve performance if hairpins are positioned after barcode sequences.}
\item{plotPositions}{logical, indicates if a density plot displaying the position of each barcode and hairpin/sgRNA sequence in the reads should be created. If \code{dualIndexForwardRead} is \code{TRUE} or \code{readfile2} is not \code{NULL}, plotPositions will generate two density plots, side by side, indicating the positions of the first barcodes and hairpins in the first plot, and second barcodes in the second.}
}

\value{
A \code{\link[edgeR:DGEList-class]{DGEList}} object with following components:
	\item{counts}{read count matrix tallying up the number of reads with particular barcode and hairpin/sgRNA matches. Each row is a hairpin/sgRNA and each column is a sample.}
	\item{genes}{In this case, hairpin/sgRNA-specific information (ID, sequences, corresponding target gene) may be recorded in this data.frame.}
	\item{lib.size}{auto-calculated column sum of the counts matrix.}
}

\details{
The \code{processAmplicons} function allows for hairpins/sgRNAs/sample index sequences to be in variable positions within each read.

The input barcode file and hairpin/sgRNA files are tab-separated text files with at least two columns (named 'ID' and 'Sequences') containing the sample or hairpin/sgRNA IDs and a second column indicating the sample index or hairpin/sgRNA sequences to be matched.
If \code{dualIndexForwardRead} is \code{TRUE}, a third column 'Sequences2' is expected in the barcode file.
If \code{readfile2} is specified, another column 'SequencesReverse' is expected in the barcode file.
The barcode file may also contain a 'group' column that indicates which experimental group a sample belongs to.
Additional columns in each file will be included in the respective \code{$samples} or \code{$genes} data.frames of the final \code{\link[edgeR:DGEList-class]{DGEList}} object.
These files, along with the FASTQ files, are assumed to be in the current working directory.

To compute the count matrix, matching to the given barcodes and hairpins/sgRNAs is conducted in two rounds.
The first round looks for an exact sequence match for the given barcode sequences and hairpin/sgRNA sequences through the entire read, returning the first match found.
If a match isn't found, the program performs a second round of matching which allows for sequence mismatches if \code{allowMismatch} is set to \code{TRUE}.
The maximum number of mismatch bases in barcode and hairpin/sgRNA are specified by the parameters \code{barcodeMismatchBase} and \code{hairpinMismatchBase} respectively.

The program outputs a \code{\link[edgeR:DGEList-class]{DGEList}} object, with a count matrix indicating the number of times each barcode and hairpin/sgRNA combination could be matched in reads from the input FASTQ files.

For further examples and data, refer to the case studies available from \url{https://bioinf.wehi.edu.au/shRNAseq/}.

The argument \code{hairpinBeforeBarcode} was introduced in edgeR release version 3.38.2 and developmental version 3.39.3.
Previously the function expected the sequences in the FASTQ files to have a fixed structure as per Figure 1A of Dai et al (2014).
The revised function can process reads where the hairpins/sgRNAs/sample index sequences are in variable positions within each read.
When \code{plotPositions=TRUE} a density plot of the  match positions is created to allow the user to assess whether they occur in the expected positions.
}

\note{
This function replaced the earlier function \code{processHairpinReads} in edgeR 3.7.17.
}

\author{Oliver Voogd, Zhiyin Dai, Shian Su and Matthew Ritchie}

\references{
Dai Z, Sheridan JM,  Gearing, LJ, Moore, DL, Su, S, Wormald, S, Wilcox, S, O'Connor, L, Dickins, RA, Blewitt, ME, Ritchie, ME (2014).
edgeR: a versatile tool for the analysis of shRNA-seq and CRISPR-Cas9 genetic screens.
\emph{F1000Research} 3, 95.
\doi{10.12688/f1000research.3928.2}
}
