getOffset <- function(y)
#	Extract offset vector or matrix from DGEList data object.
#	By default, offset is constructed from the lib.size and norm.factors
#	but offset supplied explicitly takes precedence.
#	Gordon Smyth
#	Created 26 Jan 2011. Last modified 9 Aug 2026.
{
#	Return offset component if present
	if(hasName(y,"offset")) return(y$offset)

#	Otherwise, get library sizes
	lib.size <- y$samples$lib.size
	if(is.null(lib.size)) stop("y is not a valid DGEList object")

#	Apply norm factors
	norm.factors <- y$samples$norm.factors
	if(!is.null(norm.factors)) lib.size <- lib.size*norm.factors

#	Optional prior offset defining normalization relative to the library sizes
	if(hasName(y,"offset.prior")) {
#		m <- rowMeans(y$offset.prior)
#		if(max(abs(m)) > 1e-4) y$offset.prior <- y$offset - m
		t( t(y$offset.prior) + log(lib.size) )
	} else {
		log(lib.size)
	}
}
