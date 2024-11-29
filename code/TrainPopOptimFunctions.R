# Given a certain training population size, find the set of individuals who
# 1. Maximize the lengths of marker segments IIS between the individuals and the
# two parents
# 2. Minimize the difference in those segment lengths between the two parents
# 3. Minimize the KS deviation from uniformity, done separately for each parent
# That is, the segments start at positions x. x should be uniformly distributed
# over the chromosome
# 4. Minimize the distance from the expected marker covariance matrix

# Utility assumes all individuals are inbred
trainPopUtility <- function(trainPop, parents,
                            diffLenWgt=1, ksWgt=1, covDistWgt=1,
                            chr=1, useQTL=F){
  source(here::here("code", "assessIndMrkSimilarityToPar.R"))
  source(here::here("code", "markerCovarianceDistance.R"))

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

# Not set up to go over different chromosomes, so just chr=1...
optimizeTP <- function(candidates, parents, tpSize=100,
                       diffLenWgt=1, ksWgt=1, covDistWgt=1,
                       nIter=1000,
                       useQTL=F, verbose=F){
  nCandidates <- AlphaSimR::nInd(candidates)
  tpIdx <- sample(nCandidates, size=tpSize)
  trainPop <- candidates[tpIdx]
  bestUtil <- trainPopUtility(trainPop, parents,
                              diffLenWgt=diffLenWgt, ksWgt=ksWgt, covDistWgt,
                              useQTL=useQTL) %>% sum
  allIter <- NULL
  for (optIter in 1:nIter){
    # Swap out one in trainPop
    # If new is better, keep else go back to old
    swapOut <- sample(tpSize, size=1)
    swapIn <- sample(setdiff(1:nCandidates, tpIdx), size=1)
    trainPopCand <- c(trainPop[-swapOut], candidates[swapIn])
    candUtil <- trainPopUtility(trainPopCand, parents,
                                diffLenWgt=diffLenWgt, ksWgt=ksWgt, covDistWgt,
                                useQTL=useQTL) %>% sum
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
