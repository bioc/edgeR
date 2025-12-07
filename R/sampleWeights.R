sampleWeights <- function(unit.deviance.adj, unit.df.adj, s2=NULL, iter=10L, verbose=FALSE)
# Empirical sample weights from adjusted unit deviances
# assuming that the quasi-dispersion of each observation is s2_g / w_i.
# Gordon Smyth
# Created 5 Aug 2024. Last modified 7 Dec 2025.
{
  ngenes <- nrow(unit.deviance.adj)
  nsamples <- ncol(unit.deviance.adj)
  dfgene <- rowSums(unit.df.adj)
  dfsample <- colSums(unit.df.adj)

  if(is.null(s2)) {
    w <- rep_len(1,nsamples)
    for (i in 1L:iter) {
      s2 <- drop(unit.deviance.adj %*% w) / dfgene
      s2 <- (dfgene*s2 + median(s2)) / (dfgene + 1)
      w <- colSums(unit.deviance.adj / s2) / dfsample
      w <- log(w)
      w <- exp(mean(w)-w)
      if(verbose) 
        if(nsamples > 5L)
          cat("Weights-quantiles",quantile(w,prob=c(0,0.25,0.5,0.75,1)),"\n")
        else
          cat("Weights",w,"\n")
    }
  } else {
     w <- colSums(unit.deviance.adj / s2) / dfsample
     w <- log(w)
     w <- exp(mean(w)-w)
  }

  w
}
