tpm <- function(y, effective.tx.length, rta.overdispersion=NULL, shrunk=FALSE, ...)
# Fitted tx per million from a fitted model object.
# Created 2 Oct 2025.  Last modified 2 Oct 2025.
{
  TPM <- cpm(y,log=FALSE,shrunk=shrunk,...)
  A <- effective.tx.length
  if(!is.null(rta.overdispersion)) A <- A/rta.overdispersion
  TPM <- TPM / A
  TPM / exp(mean(log(colSums(TPM)))) * 1e6
}

tpmProp <- function(tpm, geneid)
# Proportion TPM by tx for each gene
# Created 2 Oct 2025.  Last modified 2 Oct 2025.
{
  G <- rowsum(tpm,geneid,reorder=FALSE)
  tpm / G[geneid,]
}
