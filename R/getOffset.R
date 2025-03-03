getOffset <- function(y)
#	Extract offset vector or matrix from DGEList data object.
#	By default, offset is constructed from the lib.size and norm.factors
#	but offset supplied explicitly takes precedence.
#	Gordon Smyth
#	Created 26 Jan 2011. Last modified 3 Mar 2025.
{
#	Return offset if available
	if(!is.null(y$offset)) return(y$offset)

#	Otherwise, get library sizes
	lib.size <- y$samples$lib.size
	if(is.null(lib.size)) stop("y is not a valid DGEList object")

#	Apply norm factors
	norm.factors <- y$samples$norm.factors
	if(!is.null(norm.factors)) lib.size <- lib.size*norm.factors

	log(lib.size)
}
