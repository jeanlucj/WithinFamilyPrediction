# Utility function for rBayesianOptimization

# The process has a combinatorial optimization of the TP as a function of four
# parameters, though without loss of generality the first is fixed at 1
trainPopUtility <- function(diffLenWgt=0, ksWgt=0, covDistWgt=0,
                            trainPop, parents, chr=1, useQTL=F){
  # I want this to be on an exponential scale, to be able to cover more ground
  diffLenWgt <- 10^diffLenWgt
  ksWgt <- 10^ksWgt
  covDistWgt <- 10^covDistWgt

  #source(here::here("code", "assessIndMrkSimilarityToPar.R"))
  #source(here::here("code", "markerCovarianceDistance.R"))

  lenIBSp1 <- lenIBSp2 <- NULL
  for (indNum in 1:AlphaSimR::nInd(trainPop)){
    ind <- trainPop[indNum]
    lenIBSp1 <- rbind(
      lenIBSp1,
      calcMrkSimIndPar(ind, parents[1], indInbred=T, parInbred=T,
                       chr=chr, useQTL=useQTL)[[1]]
    )
    lenIBSp2 <- rbind(
      lenIBSp2,
      calcMrkSimIndPar(ind, parents[2], indInbred=T, parInbred=T,
                       chr=chr, useQTL=useQTL)[[1]]
    )
  }
  ksDev <- calcDevFromUnif(lenIBSp1) + calcDevFromUnif(lenIBSp2)
  totLen1 <- sum(lenIBSp1); totLen2 <- sum(lenIBSp2)
  diffLenPar <- abs(totLen1 - totLen2)

  # Calculate covDist, the distance between the marker covariance of the
  # training population and that expected in progeny of the F1
  empirCovMat <- calcMrkCovMat(trainPop)[[1]] # Warning: hard-coded for 1 chr
  # Warning: hard-coded for inbred parents
  cp <- matrix(1:2, nrow=1)
  f1 <- AlphaSimR::makeCross(parents, crossPlan=cp)
  exptCovMat <- calcExptMrkCovMatGamete(f1)[[1]]
  covDist <- calcMrkCovMatDist(empirCovMat, exptCovMat*4)

  return(c(totLen1 + totLen2,
           -diffLenPar * diffLenWgt,
           -ksDev * ksWgt,
           -covDist * covDistWgt))
}

makeTheTP <- function(diffLenWgt=1, ksWgt=1, covDistWgt=1,
                      candidates, parents,
                      tpSize=100, useQTL=F, chr=1,
                      nIter=1000, verbose=F){
  nCandidates <- AlphaSimR::nInd(candidates)
  tpIdx <- sample(nCandidates, size=tpSize)
  trainPop <- candidates[tpIdx]
  bestUtil <- sum(trainPopUtility(diffLenWgt=diffLenWgt,
                              ksWgt=ksWgt,
                              covDistWgt=covDistWgt,
                              trainPop, parents, chr, useQTL=useQTL))
  allIter <- NULL
  for (optIter in 1:nIter){
    # Swap out one in trainPop
    # If new is better, keep else go back to old
    swapOut <- sample(tpSize, size=1)
    swapIn <- sample(setdiff(1:nCandidates, tpIdx), size=1)
    trainPopCand <- c(trainPop[-swapOut], candidates[swapIn])
    candUtil <- sum(trainPopUtility(diffLenWgt=diffLenWgt,
                                ksWgt=ksWgt,
                                covDistWgt=covDistWgt,
                                trainPopCand, parents, chr, useQTL=useQTL))
    if (candUtil > bestUtil){
      trainPop <- trainPopCand
      tpIdx <- c(tpIdx[-swapOut], swapIn)
      bestUtil <- candUtil
    }
    allIter <- c(allIter, bestUtil)
    if (verbose) if (optIter %% 40 == 0) cat(".")
  }
  if (verbose) plot(allIter); cat("\n")
  return(trainPop)
}
