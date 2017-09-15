cpmByGroup <- function(y, ...)
UseMethod("cpmByGroup")

cpmByGroup.DGEList <- function(y, group=NULL, log=FALSE, prior.count=0.25, ...)
#	Counts per million averaged by group
#	Gordon Smyth
#	Created 10 July 2017
{
	if(!is.null(group)) y$samples$group <- group
	if(!log) prior.count <- 0
	fit <- glmFit(y,prior.count=prior.count,...)

	logCPM <- fit$coefficients + log(1e6)
	if(log) {
		logCPM/log(2)
	} else {
		exp(logCPM)
	}
}

cpmByGroup.default <- function(y, group=NULL, dispersion=0.05, log=FALSE, prior.count=0.25, ...)
#	Counts per million averaged by group
#	Gordon Smyth
#	Created 10 July 2017
{
	if(!log) prior.count <- 0
	if(is.null(group)) {
		fit <- glmFit(y,dispersion=dispersion,prior.count=prior.count,...)
	} else {
		group <- as.factor(group)
		design <- model.matrix(~0+group)
		colnames(design) <- levels(group)
		fit <- glmFit(y,design=design,dispersion=dispersion,prior.count=prior.count,...)
	}

	logCPM <- fit$coefficients + log(1e6)
	if(log) {
		logCPM/log(2)
	} else {
		exp(logCPM)
	}
}
