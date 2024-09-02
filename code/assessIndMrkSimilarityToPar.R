# Count the number of sequential markers that have the same states as the parent
# What the function needs
# 1. The parent, so it can retrieve the parent haplotypes
# 2. Whether the parent is completely inbred or not
# 3. What chromosome it's working on
# 4. The individual
# 5. Whether the individual is completely inbred or not
# I'm going to pass both the individual and the parent as Formal class pop
# Return a list with one, two or four vectors
# Each vector contains the number of markers in the sequence, its start position
# As of now, nSnpPerChr will have to be a global variable
calcMrkSimIndPar <- function(ind, par, indInbred=T, parInbred=T, chr=1, useQTL=F){
  if (useQTL){
    indMrk <- AlphaSimR::pullSegSiteHaplo(ind, chr=chr)
    parMrk <- AlphaSimR::pullSegSiteHaplo(par, chr=chr)
  } else{
    mrkNames <- paste(chr, 1:nSnpPerChr, sep="_")
    indMrk <- AlphaSimR::pullMarkerHaplo(ind, markers=mrkNames)
    parMrk <- AlphaSimR::pullMarkerHaplo(par, markers=mrkNames)
  }
  toRet <- list()
  idx <- 1
  for (indHaplo in 1:(2 - indInbred)){
    for (parHaplo in 1:(2 - parInbred)){
      toRet[[idx]] <- findPosMaxSeqSame(indMrk[indHaplo,], parMrk[parHaplo,])
      idx <- idx + 1
    }
  }
  return(toRet)
}

# Find the longest number of contiguous markers that are IBS
# If there are 2 or more longest sequences of the same length, this only reports
#  the first. I need to fix this to enable it to report _all_ max length seq.
# Return two integers: the number of markers in the sequence, its start position
findPosMaxSeqSame <- function(parHaplo, indHaplo){
  nMrk <- length(parHaplo)
  if (nMrk != length(indHaplo))
    stop("findPosMaxSeqSame: parHaplo and indHaplo should be same length")
  mrkSame <- parHaplo == indHaplo
  maxSame <- 0
  strtPos <- 0
  while (strtPos < nMrk){
    lengthSame <- 0
    pos <- strtPos
    while (mrkSame[pos+1]){
      lengthSame <- lengthSame+1
      pos <- pos+1
      if (lengthSame > maxSame){
        maxSame <- lengthSame
        posMaxSame <- strtPos+1
      }
      if (pos == nMrk) break
    }
    strtPos <- pos+1
  }
  return(c(maxSame, posMaxSame))
}

# Propose a training population
# Each member of this population will have a longest matching sequence that
# starts at position x. Assemble into a vector and calculate the deviation of
# the empirical distribution from uniform (this could go into a Kolmogorov -
# Smirnov test)
# matchSeqMat has the length of matches in col 1 and the start position in
# col 2
# nMrk is the number of markers: the uniform goes from 1 to nMrk
calcDevFromUnif <- function(matchSeqMat, nMrk=nSnpPerChr){
  # pool together start and end position of sequences (uniform together)
  nInd <- nrow(matchSeqMat)
  strtEnd <- c(matchSeqMat[,2], matchSeqMat[,2] + matchSeqMat[,1] - 1) %>% sort
  expected <- 1 + (nMrk - 1) * (2 * 1:(nInd*2) - 1) / (nInd*4)
  return(max(strtEnd - expected))
}

# Hmmm. I don't actually need this.
# This is some version of BLAST (I think). I'm starting with an initial hash of
# three markers sequences
# haplo is a vector of 0 and 1 allelic states for markers
makeMrkHash <- function(haplo){
  nMrk <- length(haplo)
  hash <- list(integer(0))
  for (i in 1:7) hash <- c(hash, list(integer(0)))
  for (strt in 1:(nMrk - 2)){
    mrkSeq <- haplo[strt + 0:2]
    idx <- 1 + mrkSeq %*% c(4, 2, 1)
    hash[[idx]] <- c(hash[[idx]], strt)
  }
  return(hash)
}
