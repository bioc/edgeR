PCList <- function(counts, counts2, samples=NULL, group=NULL, genes=NULL, ...)
# Created by Lizhong Chen
#	11 April 2025.  Last modified 15 April 2025.
{
#	Check counts
	counts  <- as.matrix(counts)
	counts2 <- as.matrix(counts2)
	
	if(!length(counts)) stop("'counts' must contain at least one value")
	if(!length(counts2)) stop("'counts2' must contain at least one value")

	if(!identical(dim(counts),dim(counts2))) stop("The dimensions of 'counts' and 'counts2' are not matched")
	
	m1 <- min(counts)
	m2 <- min(counts2)
	if(is.na(m1) || is.na(m2)) stop("NA counts not allowed")
	if(m1 < 0 || m2 < 0) stop("Negative counts not allowed")
	if(is.infinite(max(counts)) || is.infinite(max(counts2))) stop("Infinite counts not allowed")
	
	nlib  <- ncol(counts)
	ntags <- nrow(counts)
	
	if(is.null(colnames(counts))) colnames(counts) <- colnames(counts2) <- paste0("Sample",1:nlib)
	if(is.null(rownames(counts))) rownames(counts) <- rownames(counts2) <- 1:ntags

#	Check samples
	if(!is.null(samples)) {
		samples <- as.data.frame(samples)
		if(nlib != nrow(samples)) stop("Number of rows in 'samples' must equal number of columns in 'counts'")
	}
	
#	Get group from samples if appropriate
	if(is.null(group) && !is.null(samples$group)) {
		group <- samples$group
		samples$group <- NULL
	}
	
#	Check group
	if(is.null(group)) {
		group <- rep_len(1L,nlib)
		levels(group) <- "1"
		class(group)  <- "factor"
	} else {
		if(length(group) != nlib) stop("Length of 'group' must equal number of columns in 'counts'")
		group <- dropEmptyLevels(group)
	}
	
#	Make data frame of sample information
	sam <- data.frame(group=group)
	if(!is.null(samples)) sam <- data.frame(sam, samples)
	samples <- sam
	if(anyDuplicated(colnames(counts))) {
		message("Repeated column names found in count matrix")
		row.names(samples) <- 1L:nlib
	} else 
		row.names(samples) <- colnames(counts)
	
#	Make object
	x <- new("PCList", list(counts=counts, counts2=counts2, samples=samples))
	
#	Add data frame of gene information
	if(!is.null(genes)) {
		genes <- as.data.frame(genes, stringsAsFactors=FALSE)
		if(nrow(genes) != ntags) stop("Counts and genes have different numbers of rows")
		if(anyDuplicated(row.names(counts)))
			warning("Count matrix has duplicated rownames",call.=FALSE)
		else 
			row.names(genes) <- row.names(counts)
		x$genes <- genes
	}
	
	x
}