# Haldane's mapping function relating the distance in Morgans d to the
# recombination rate r
# r = (1 / 2) (1 - e^(-2d))
# d = -(1 / 2) log(1 - 2r)
# The covariance is +/- 0.5*(1 - r) - 0.25 = 0.25 - 0.5*r
# AlphaSimR uses some other mapping function function McPeek and Speed 1995
# But I'm just going to assume that the marker positions are distances in
# Morgans and use Haldane's mapping function because that's easy
#
# Function to calculate the covariance matrix across markers in a population
calcMrkCovMat <- function(pop){
  snpMap <- AlphaSimR::getSnpMap()
  snps <- AlphaSimR::pullSnpGeno(pop)
  allChrCov <- list()
  for (chr in unique(snpMap$chr)){
    snpsOnChr <- snps[,snpMap$chr == chr]
    allChrCov <- c(allChrCov, list(cov(snpsOnChr)))
  }
  return(allChrCov)
}

# Function to calculate the expected covariance among markers on gametes from
# an individual
calcExptMrkCovMatGamete <- function(ind){
  snpMap <- AlphaSimR::getSnpMap()
  snps <- AlphaSimR::pullSnpHaplo(ind)*2L - 1L
  allChrCov <- list()
  for (chr in unique(snpMap$chr)){
    snpsOnChrPos <- snpMap$pos[snpMap$chr == chr]
    distMat <- dist(snpsOnChrPos, diag=T, upper=T)
    recFreqMat <- as.matrix(0.5 * (1 - exp(-2*distMat)))
    rOnCov <- 0.25 - 0.5*recFreqMat
    snpsOnChr <- snps[,snpMap$chr == chr]
    # If exptCovMat row col is + the covariance is + else the covariance is -
    exptCovMat <- tcrossprod(snpsOnChr[1,]) * rOnCov
    monomorphic <- snpsOnChr[1,] == snpsOnChr[2,]
    exptCovMat[monomorphic,] <- exptCovMat[,monomorphic] <- 0
    allChrCov <- c(allChrCov, list(exptCovMat))
  }
  return(allChrCov)
}

# Calculates the Manhattan distance of the upper.tri including the diagonal
calcMrkCovMatDist <- function(covMat1, covMat2){
  require(magrittr)
  if (any(dim(covMat1) != dim(covMat2)))
    stop("calcMrkCovMatDist: Matrices have to be the same size")
  if (length(unique(dim(covMat1))) != 1)
    stop("calcMrkCovMatDist: Matrices have to be square")
  coef1 <- covMat1[upper.tri(covMat1, diag=T)]
  coef2 <- covMat2[upper.tri(covMat2, diag=T)]
  return(sum(abs(coef1 - coef2)))
}

