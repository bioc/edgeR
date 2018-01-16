nearestTSS <- function(chr,locus,species="Hs")
#	Find nearest gene transcriptional start sites from orgDb
#	Gordon Smyth
#	Created 3 Jan 2018.  Last modified 11 Jan 2018.
{
#	Get access to required annotation functions
	suppressPackageStartupMessages(OK <- requireNamespace("AnnotationDbi",quietly=TRUE))
	if(!OK) stop("AnnotationDbi package required but not installed (or can't be loaded)")

#	Load appropriate organism package
	orgPkg <- paste0("org.",species,".eg.db")
	suppressPackageStartupMessages(OK <- requireNamespace(orgPkg,quietly=TRUE))
	if(!OK) stop(orgPkg," package required but not not installed (or can't be loaded)")

#	Get gene start positions
	obj <- paste0("org.",species,".egCHRLOC")
	egCHRLOC <- tryCatch(getFromNamespace(obj,orgPkg), error=function(e) FALSE)
	if(is.logical(egCHRLOC)) stop("Can't find egCHRLOC gene location mappings in package ",orgPkg)
	EGLOC <- AnnotationDbi::toTable(egCHRLOC)

#	Get first gene end position
	obj <- paste0("org.",species,".egCHRLOCEND")
	egCHRLOCEND <- tryCatch(getFromNamespace(obj,orgPkg), error=function(e) FALSE)
	if(is.logical(egCHRLOCEND)) stop("Can't find egCHRLOCEND gene end mappings in package ",orgPkg)
	EGEND <- AnnotationDbi::toTable(egCHRLOCEND)
	EGLOC$end_location <- EGEND$end_location

#	Get Symbols
	obj <- paste0("org.",species,".egSYMBOL")
	egSYMBOL <- tryCatch(getFromNamespace(obj,orgPkg), error=function(e) FALSE)
	if(is.logical(egSYMBOL)) stop("Can't find egSYMBOL gene symbol mappings in package ",orgPkg)
	EGSym <- AnnotationDbi::toTable(egSYMBOL)
	m <- match(EGLOC$gene_id,EGSym$gene_id)
	EGLOC$symbol <- EGSym[m,2]

#	Get strand, width and TSS
	EGLOC$neg <- (EGLOC$start_location < 0L)
	EGLOC$width <- EGLOC$end_location - EGLOC$start_location
	EGLOC$tss <- EGLOC$start_location+1L
	EGLOC$tss[EGLOC$neg] <- EGLOC$end_location[EGLOC$neg]

#	Keep only positive loci
	EGLOC$width <- abs(EGLOC$width)
	EGLOC$tss <- abs(EGLOC$tss)
	EGLOC$start_location <- EGLOC$end_location <- NULL
	EGLOC$strand <- rep_len("+",nrow(EGLOC))
	EGLOC$strand[EGLOC$neg] <- "-"

#	Sort by genomic position
	o <- order(EGLOC$Chromosome,EGLOC$tss)
	EGLOC <- EGLOC[o,]
#	NREF <- nrow(EGLOC)
#	ChrFirst <- which(EGLOC$Chromosome[-NREF] != EGLOC$Chromosome[-1L])
#	ChrN <- c(1L,ChrFirst+1L)

#	Do chr values start with "chr"?
	if(length(grep("^chr",chr[1]))) EGLOC$Chromosome <- paste0("chr",EGLOC$Chromosome)

#	Prepare output
	ILocus <- rep_len(0L,length(chr))

#	Cycle over chromosomes
	ChrNames <- unique(EGLOC$Chromosome)
	for (ChrA in ChrNames) {
		iref <- which(EGLOC$Chromosome==ChrA)
		iinc <- which(chr==ChrA)
		Which <- nearestReftoX(locus[iinc], EGLOC$tss[iref])
		ILocus[iinc] <- iref[Which]
	}

#	Expand EGLOC to input length
	EGLOC$Chromosome <- NULL
	ILocus[ILocus==0L] <- NA
	Out <- EGLOC[ILocus,,drop=FALSE]
	Out$distance <- Out$tss - locus
	Out$distance[Out$neg] <- -Out$distance[Out$neg]
	Out$neg <- NULL

	Out
}
