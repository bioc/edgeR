catchRSEM <- function(paths,nboot=100,verbose=TRUE)
  # Read transcriptwise counts and Gibbs samples from RSEM output
  # Use Gibbs samples to estimate overdispersion of transcriptwise counts
  # Pedro Baldoni and Gordon Smyth
  # Created 24 April 2024. Last modified 11 May 2024.
{
  NSamples <- length(paths)
  NBoot <- nboot
  
  OK <- requireNamespace("readr",quietly=TRUE)
  if(!OK) stop("readr package required but is not installed (or can't be loaded)")
  
  # Initialize vector of inferential sample types
  Type <- "gibbs"
  ResampleType <- rep_len(Type,NSamples)
  
  # Accumulate counts and CV^2 of bootstrap counts for each sample
  for (j in 1L:NSamples) {
    if(verbose) cat("Reading ",paths[j],", ",sep="")
    
    #  File locations
    QuantFile <- list.files(paths[j],"isoforms.results",full.names = TRUE)
    if(!file.exists(QuantFile)) stop("isoforms.results file not found at specified path", call.=FALSE)
    if(length(QuantFile)>1) stop("More than one isoforms.results file found at specified path", call.=FALSE)
    
    #  Read counts
    if(j == 1L) {
      Quant1 <- suppressWarnings(readr::read_tsv(QuantFile,col_types="c_ddd___dd___",progress=FALSE))
      NTx <- nrow(Quant1)
      Counts <- matrix(0,NTx,NSamples)
      DF <- rep_len(0L,NTx)
      OverDisp <- rep_len(0,NTx)
      Counts[,1L] <- Quant1$expected_count
      M <- Quant1$posterior_mean_count
      S <- Quant1$posterior_standard_deviation_of_count
    } else {
      Quant <- suppressWarnings(readr::read_tsv(QuantFile,col_types="____d___dd___",progress=FALSE))
      Counts[,j] <- Quant$expected_count
      M <- Quant$posterior_mean_count
      S <- Quant$posterior_standard_deviation_of_count
    }
    
    if(verbose) cat(NTx,"transcripts,",NBoot,Type,"samples\n")
    
    #  Bootstrap samples
    if(NBoot > 0L) {
      i <- (M > 0)
      OverDisp[i] <- OverDisp[i] + (NBoot-1)*(S[i]^2)/M[i]
      DF[i] <- DF[i]+NBoot-1L
    }
  }
  
  # Estimate overdispersion for each transcript
  i <- (DF > 0L)
  if(sum(i) > 0L) {
    OverDisp[i] <- OverDisp[i] / DF[i]
    #  Apply a limited amount of moderation
    DFMedian <- median(DF[i])
    DFPrior <- 3
    OverDispPrior <- median(OverDisp[i]) / qf(0.5,df1=DFMedian,df2=DFPrior)
    if(OverDispPrior < 1) OverDispPrior <- 1
    OverDisp[i] <- (DFPrior * OverDispPrior + DF[i]*OverDisp[i]) / (DFPrior + DF[i])
    OverDisp <- pmax(OverDisp,1)
    OverDisp[!i] <- OverDispPrior
  } else {
    OverDisp[] <- NA_real_
    OverDispPrior <- NA_real_
  }
  
  # Prepare output
  Quant1 <- as.data.frame(Quant1,stringsAsFactors=FALSE)
  dimnames(Counts) <- list(Quant1$transcript_id,paths)
  row.names(Quant1) <- Quant1$transcript_id
  Quant1$transcript_id <- Quant1$expected_count <- NULL
  Quant1$posterior_mean_count <- Quant1$posterior_standard_deviation_of_count<- NULL
  colnames(Quant1) <- c('Length','EffectiveLength')
  Quant1$Overdispersion <- OverDisp
  
  list(counts=Counts,annotation=Quant1,overdispersion.prior=OverDispPrior,resample.type=ResampleType)
}
