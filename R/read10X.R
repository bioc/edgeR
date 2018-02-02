read10X <- function(mtx="matrix.mtx",genes="genes.tsv",barcodes=NULL,path=NULL,DGEList=TRUE)
{
#	Add path
	if(!is.null(path)) {
		mtx <- file.path(path,mtx)
		genes <- file.path(path,genes)
		if(!is.null(barcodes)) barcodes <- file.path(path,barcodes)
	}

#	Fetch header info for checking
	N <- scan(mtx,skip=2,what=0L,sep=" ",n=3,quiet=TRUE)
	ngenes <- N[1]
	ncells <- N[2]
	nmtx <- N[3]

#	Read gene Ids
	Genes <- read.table(genes,header=FALSE,comment.char="",sep="\t",row.names=1,colClasses="character")
	names(Genes) <- "Symbol"
	if(nrow(Genes) != ngenes) stop("Number of gene IDs doesn't agree with header information in mtx file")

#	Read mtx file of counts
	m <- read.table(mtx,skip=3,header=FALSE,comment.char="",sep=" ",colClasses="integer",nrows=nmtx)

#	Convert Market Exchange Format to ordinary matrix
	y <- matrix(0L,ngenes,ncells)
	i <- m[,1]+(m[,2]-1L)*ngenes
	y[i] <- m[,3]
	dimnames(y) <- list(Gene=row.names(Genes),Cell=1:ncells)

#	Optionally read barcodes
	if(is.null(barcodes)) {
		Samples <- NULL
	} else {
		Barcodes <- scan(barcodes,what="",quiet=TRUE)
		if(length(Barcodes) != ncells) stop("Number of barcodes doesn't agree with header information in mtx file")
		Samples <- data.frame(Barcode=Barcodes)
	}

	if(DGEList) {
		DGEList(count=y,genes=Genes,samples=Samples)
	} else {
		list(counts=y,samples=Samples,genes=Genes)
	}
}
