normalizeBetweenArrays.DGEList <- function(object, method = "cyclicloess", cyclic.method = "affy", ...)
# Apply limma normalization methods to DGEList object via the offset matrix.
# Gordon Smyth.
# Created 10 Aug 2024. Last modified 11 Aug 2024.
{
  method <- match.arg(method, c("none","scale","quantile","cyclicloess"))
  logCPM.raw <- cpm(object, log=TRUE)
  logCPM.norm <- normalizeBetweenArrays(logCPM.raw, method=method, cyclic.method=cyclic.method, ...)
  object$offset <- matrix(log(object$samples$lib.size), nrow(object), ncol(object), byrow=TRUE)
  object$offset <- object$offset + (logCPM.raw - logCPM.norm) * log(2)
  object
}
