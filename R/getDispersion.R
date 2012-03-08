getDispersion <- function(y)
#	Get most complex dispersion values from DGEList object
#	Gordon Smyth
#	Created 12 Dec 2011.  Last modified 1 March 2012.
{
	if( !is.null(y$tagwise.dispersion) ) {
		dispersion <- y$tagwise.dispersion
		type <- "tagwise"
	} else {
		if( !is.null(y$trended.dispersion) ) {
			dispersion <- y$trended.dispersion
			type <- "trended"
		} else {
			if( !is.null(y$common.dispersion) ) {
				dispersion <- y$common.dispersion
				type <- "common"
			} else
				dispersion <- NULL
		}
	}
	dispersion
}