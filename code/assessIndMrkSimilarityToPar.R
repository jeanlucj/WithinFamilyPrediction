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
calcMrkSimIndPar <- function(ind, par, indInbred=T, parInbred=T,
                             chr=1, useQTL=F){
  if (useQTL){
    indMrk <- AlphaSimR::pullSegSiteHaplo(ind, chr=chr)
    parMrk <- AlphaSimR::pullSegSiteHaplo(par, chr=chr)
  } else{
    indMrk <- AlphaSimR::pullSnpHaplo(ind, chr=chr)
    parMrk <- AlphaSimR::pullSnpHaplo(par, chr=chr)
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
      if (pos == nMrk) break
    }
    if (lengthSame > maxSame){
      maxSame <- lengthSame
      posMaxSame <- strtPos+1
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
# nSegSitePerChr is the number of markers: the uniform goes from 1 to
# nSegSitePerChr
# NOTE:
# If you are using the QTL, nSegSitePerChr = nQTL + nSNP
# If you are NOT using the QTL, nSegSitePerChr = nSNP
calcDevFromUnif <- function(matchSeqMat, nSegSitePerChr=200){
  # pool together start and end position of sequences (uniform together)
  nInd <- nrow(matchSeqMat)
  strtEnd <- sort(c(matchSeqMat[,2], matchSeqMat[,2] + matchSeqMat[,1] - 1))
  expected <- 1 + (nSegSitePerChr - 1) * (2 * 1:(nInd*2) - 1) / (nInd*4)
  return(max(strtEnd - expected))
}

# Moving from left to right along a chromosome, if two haplotypes are the same
# at the first locus, what probability are they the same at the second locus?
# Calculate those probabilities for all loci and all chromosomes Returns a list
# of three-row matrices.  The list is nChr long. The first row of each matrix is
# the frequency of the 1 allele, the second row is if the allele at the first
# locus is 0, the second row is if the allele is 1.
calcProbSameSame <- function(pop, useQTL=F){
  nInd <- AlphaSimR::nInd(pop)
  nChr <- pop@nChr
  probSameSame <- list()
  for (chr in 1:nChr){
    if (useQTL){
      mrkDat <- AlphaSimR::pullSegSiteHaplo(pop, chr=chr)
    } else{
      mrkDat <- AlphaSimR::pullSnpHaplo(pop, chr=chr)
    }
    nMrkPerChr <- ncol(mrkDat)
    # Calculates both the allele frequency of the "1" allele and the probability
    # that the next locus is IIS given that this locus is and the allele is 0
    # and the next locus is IIS given that this locus is and the allele is 1
    # So, three numbers...
    calcPSSoneMrk <- function(mrk){
      probAllele1 <- mean(mrkDat[,mrk])
      rstuTab <- table(2*mrkDat[,mrk] + mrkDat[,mrk+1]) / nInd / 2
      rstu <- numeric(4)
      names(rstu) <- as.character(0:3)
      rstu[names(rstuTab)] <- rstuTab
      return(c(probAllele1,
               sum(rstu[1:2]^2) / sum(rstu[1:2])^2,
               sum(rstu[3:4]^2) / sum(rstu[3:4])^2))
    }
    probSameSame <- c(probSameSame,
                      list(
                        cbind(
                          sapply(1:(nMrkPerChr-1), calcPSSoneMrk),
                          c(mean(mrkDat[,nMrkPerChr]), NA, NA)
                          )))
  }
  return(probSameSame)
}

# This function is analogous to findPosMaxSeqSame except that instead of
# counting the number of contiguous markers that are IIS it calculates the
# probability of the IIS sequence and returns the position of the sequence with
# the lowest probability, the thought being that low probability would not be
# assembled by chance but would have arisen by transmission from a recent common
# ancestor
# Find the sequence of contiguous markers that are IIS with the lowest probability
# Return two values: -log10 of the probability, the sequence start position
findPosMinProbSame <- function(parHaplo, indHaplo, probSmSmChr){
  nMrk <- length(parHaplo)
  if (nMrk != length(indHaplo))
    stop("findPosMaxSeqSame: parHaplo and indHaplo should be same length")
  mrkSame <- parHaplo == indHaplo
  minProbSame <- 1
  strtPos <- 0
  while (strtPos < nMrk){
    pos <- strtPos
    logProbSame <- 0
    while (mrkSame[pos+1]){
      if (logProbSame == 0){
        alleleFreq <- probSmSmChr[1, pos+1]
        if (indHaplo[pos+1] == 0) alleleFreq <- 1 - alleleFreq
        logProbSame <- log(alleleFreq)
      } else{
        logProbSame <- logProbSame + log(probSmSmChr[2+indHaplo[pos], pos])
      }
      pos <- pos+1
      if (pos == nMrk) break
    }
    if (logProbSame < minProbSame){
      minProbSame <- logProbSame
      posMinProbSame <- strtPos+1
    }
    strtPos <- pos+1
  }
  return(c(minProbSame, posMinProbSame))
}

# This function calculates the probabilities of ALL streaks of contiguous
# markers that are IIS. It returns the position of the streak, its length, and
# its probability for all such streaks on a _chromosome_. Return matrix with
# three columns: -log of streak probability, position of the streak, its length
# Note: this also works if you want to use the QTL. Then calculate probSmSmChr
# with useQTL=T and include QTL in parHaplo and indHaplo.
calcAllPosProbSame <- function(parHaplo, indHaplo, probSmSmChr){
  nMrk <- length(parHaplo)
  if (nMrk != length(indHaplo))
    stop("calcAllPosProbSame: parHaplo and indHaplo should be same length")
  mrkSame <- parHaplo == indHaplo
  toRet <- NULL
  strtPos <- 0
  while (strtPos < nMrk){
    pos <- strtPos
    logProbSame <- 0
    while (mrkSame[pos+1]){
      if (logProbSame == 0){
        alleleFreq <- probSmSmChr[1, pos+1]
        if (indHaplo[pos+1] == 0) alleleFreq <- 1 - alleleFreq
        logProbSame <- log(alleleFreq)
      } else{
        logProbSame <- logProbSame + log(probSmSmChr[2+indHaplo[pos], pos])
      }
      pos <- pos+1
      if (pos == nMrk) break
    }
    if (logProbSame != 0){
      toRet <- rbind(toRet, c(logProbSame, strtPos+1, pos - strtPos))
    }
    strtPos <- pos+1
  }
  return(toRet)
}

# Function to split markers into one matrix if they appear to be IBD with either
# parent of a biparental and a separate matrix if they they are not.
# pop is the overall population that has both genotypes and phenotypes
# parents is a two vector with indices of the parents that are assumed in pop
# threshold is a scalar with the probability required to put into the one
# matrix versus the other, one threshold for each parent
# For now (29 Sept 2024) I'm going to assume all individuals are fully inbred
partitionMrkForBiparPred <- function(pop, parents, threshold){
  pss <- calcProbSameSame(pop)
  snpDat <- AlphaSimR::pullSnpHaplo(pop)
  nInd <- AlphaSimR::nInd(pop)
  snpMap <- AlphaSimR::getSnpMap()
  nChr <- pop@nChr
  nSnp <- nrow(snpMap)
  ibdMat <- nonIbdMat <- matrix(nrow=nInd, ncol=nSnp)
  for (ind in 1:nInd){
    for (chr in 1:nChr){
      chrSnp <- which(snpMap$chr == chr)
      toIbdMat <- NULL
      for (par in 1:2){
        parIdx <- parents[par]
        # The *2-1 calculation is to skip every other row since inbred
        chrProbs <- calcAllPosProbSame(snpDat[parIdx*2-1, chrSnp],
                                       snpDat[ind*2-1, chrSnp], pss[[chr]])
        chrProbs <- chrProbs[chrProbs[,1] < threshold,,drop=F]
        if (nrow(chrProbs) > 0){
          for (segment in 1:nrow(chrProbs)){
            toIbdMat <- c(toIbdMat,
                          chrProbs[segment, 2] - 1 + 1:chrProbs[segment, 3])
          }
        }
      }
      if (!is.null(toIbdMat)){
        toIbdMat <- chrSnp[toIbdMat]
        ibdMat[ind, toIbdMat] <- snpDat[ind*2-1, toIbdMat]
      }
      toNonIbdMat <- setdiff(chrSnp, toIbdMat)
      nonIbdMat[ind, toNonIbdMat] <- snpDat[ind*2-1, toNonIbdMat]
    }#END chr
  }#END ind
  # Fill in the NA of these matrices
  fillNAwithMean <- function(v){
    v[is.na(v)] <- mean(v, na.rm=T)
    return(v)
  }
  nonIbdMat <- apply(nonIbdMat, 2, fillNAwithMean)
  ibdMat <- apply(ibdMat, 2, fillNAwithMean)
  return(list(ibdMat=ibdMat, nonIbdMat=nonIbdMat))
}
