DGEList <- function(counts, ...)
UseMethod("DGEList")

DGEList.default <- function(counts, lib.size=NULL, norm.factors=NULL, samples=NULL, group=NULL, genes=NULL, remove.zeros=FALSE, ...)
#	Construct DGEList object from components, with some checking
#	Created 28 Sep 2008. Last modified 1 Mar 2023.
{
#	Check counts
	counts <- as.matrix(counts)
	if(!length(counts)) stop("'counts' must contain at least one value")
	m <- min(counts)
	if(is.na(m)) stop("NA counts not allowed")
	if(m < 0) stop("Negative counts not allowed")
	if(is.infinite(max(counts))) stop("Infinite counts not allowed")
	nlib <- ncol(counts)
	ntags <- nrow(counts)
	if(is.null(colnames(counts))) colnames(counts) <- paste0("Sample",1:nlib)
	if(is.null(rownames(counts))) rownames(counts) <- 1:ntags

#	Check lib.size
	if(is.null(lib.size)) {
		lib.size <- colSums(counts)
		if(min(lib.size) <= 0) warning("At least one library size is zero")
	} else {
		if(!is.numeric(lib.size)) stop("'lib.size' must be numeric")
		if(!identical(length(lib.size),nlib)) stop("length of 'lib.size' must equal number of samples")
		m <- min(lib.size)
		if(is.na(m)) stop("NA library sizes not allowed")
		if(m < 0) stop("negative library sizes not allowed")
		if(identical(m,0)) {
			if(any(lib.size==0 & colSums(counts)>0)) stop("library size set to zero but counts for that sample are nonzero")
		}
	}

#	Check norm.factors
	if(is.null(norm.factors)) {
		norm.factors <- rep_len(1,ncol(counts))
	} else {
		if(!is.numeric(norm.factors)) stop("'lib.size' must be numeric")
		if(!identical(nlib,length(norm.factors))) stop("Length of 'norm.factors' must equal number of columns in 'counts'")
		m <- min(norm.factors)
		if(is.na(m)) stop("NA norm factors not allowed")
		if(m <= 0) stop("norm factors must be positive")
		if( abs(sum(log(norm.factors)) > 1e-6 )) warning("norm factors don't multiply to 1")
	}

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
		class(group) <- "factor"
	} else {
		if(length(group) != nlib) stop("Length of 'group' must equal number of columns in 'counts'")
		group <- dropEmptyLevels(group)
	}

#	Make data frame of sample information
	sam <- data.frame(group=group,lib.size=lib.size,norm.factors=norm.factors)
	if(!is.null(samples)) sam <- data.frame(sam, samples)
	samples <- sam
	if(anyDuplicated(colnames(counts))) {
		message("Repeated column names found in count matrix")
		row.names(samples) <- 1L:nlib
	} else 
		row.names(samples) <- colnames(counts)

#	Make object
	x <- new("DGEList",list(counts=counts,samples=samples))

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

#	Remove rows with all zeros
	if(remove.zeros) {
		all.zeros <- rowSums(counts>0,na.rm=TRUE)==0
		if(any(all.zeros)) {
			x <- x[!all.zeros,]
			message("Removing ",sum(all.zeros)," rows with all zero counts")
		}
	}

#	x$offset <- expandAsMatrix(getOffset(x),dim(counts))
#	x$weights <- matrix(1,ntags,nlib)

	x
}

DGEList.data.frame <- function(counts, lib.size=NULL, norm.factors=NULL, samples=NULL, group=NULL, genes=NULL, remove.zeros=FALSE, annotation.columns=NULL, ...)
#	Construct DGEList object from components, with some checking
#	Gordon Smyth
#	Created 28 Feb 2023. Last modified 1 Mar 2023.
{
	if(is.null(annotation.columns)) {
#		Annotation columns are unspecified, but check for undeclared annotation columns
		cl <- vapply(counts,class,FUN.VALUE="")
		IsNumeric <- (cl %in% c("integer","numeric"))
		NColumns <- length(IsNumeric)
		AlphaColumns <- which(!IsNumeric)
		if(length(AlphaColumns)) {
			LastAlphaColumn <- AlphaColumns[length(AlphaColumns)]
			if(identical(LastAlphaColumn,NColumns)) stop("`counts` appears to have no numeric columns.")
			annotation.columns <- 1:LastAlphaColumn
			numeric.columns <- (LastAlphaColumn+1L):NColumns
			if(identical(LastAlphaColumn,1L))
				message("Setting first column of `counts` as gene annotation.")
			else
				message("Non-numeric columns found: setting first ",LastAlphaColumn," columns of `counts` to be gene annotation.")
			if(is.null(genes))
				genes <- counts[,annotation.columns,drop=FALSE]
			else
				genes <- data.frame(counts[,annotation.columns,drop=FALSE],genes)
			counts <- counts[,numeric.columns,drop=FALSE]
		}
	} else {
#		If annotation columns are specified, separate them out from the numeric columns of counts
#		First convert annotation.columns to integers
		if(is.logical(annotation.columns)) annotation.columns <- which(annotation.columns)
		if(is.character(annotation.columns)) annotation.columns <- which(names(counts) %in% annotation.columns)
		if(is.null(genes))
			genes <- counts[,annotation.columns,drop=FALSE]
		else
			genes <- data.frame(counts[,annotation.columns,drop=FALSE],genes)
		counts <- counts[,-annotation.columns,drop=FALSE]
	}

#	Call default method
	counts <- as.matrix(counts)
	DGEList(counts=counts,lib.size=lib.size,norm.factors=norm.factors,samples=samples,group=group,genes=genes,remove.zeros=remove.zeros)
}

.isAllZero <- function(y) 
# Check whether all counts are zero.
# Also checks and stops with an informative error message if negative, NA or infinite counts are present.
{
	if (!length(y)) return(FALSE)
	check.range <- range(y)
	if (is.na(check.range[1])) stop("NA counts not allowed", call.=FALSE)
	if (check.range[1] < 0) stop("Negative counts not allowed", call.=FALSE)
	if (is.infinite(check.range[2])) stop("Infinite counts not allowed", call.=FALSE)
	check.range[2]==0
}
