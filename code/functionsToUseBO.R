# Utility function for rBayesianOptimization

# The process has a combinatorial optimization of the TP as a function of four
# parameters, though without loss of generality the first is fixed at 1
trainPopUtility <- function(diffLenWgt=0, ksWgt=0, covDistWgt=0,
                            trainPop, parents, chr=1, useQTL=F){
  # I want this to be on an exponential scale, to be able to cover more ground
  diffLenWgt <- 10^diffLenWgt
  ksWgt <- 10^ksWgt
  covDistWgt <- 10^covDistWgt

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

makeTheTP <- function(diffLenWgt=1, ksWgt=1, covDistWgt=1,
                       candidates, parents,
                       tpSize=100, useQTL=F, chr=1,
                       nIter=1000, verbose=F){
  nCandidates <- AlphaSimR::nInd(candidates)
  tpIdx <- sample(nCandidates, size=tpSize)
  trainPop <- candidates[tpIdx]
  bestUtil <- trainPopUtility(diffLenWgt=diffLenWgt,
                              ksWgt=ksWgt,
                              covDistWgt=covDistWgt,
                              trainPop, parents, chr, useQTL=useQTL) %>% sum
  allIter <- NULL
  for (optIter in 1:nIter){
    # Swap out one in trainPop
    # If new is better, keep else go back to old
    swapOut <- sample(tpSize, size=1)
    swapIn <- sample(setdiff(1:nCandidates, tpIdx), size=1)
    trainPopCand <- c(trainPop[-swapOut], candidates[swapIn])
    candUtil <- trainPopUtility(diffLenWgt=diffLenWgt,
                                ksWgt=ksWgt,
                                covDistWgt=covDistWgt,
                                trainPopCand, parents, chr, useQTL=useQTL) %>% sum
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

# I don't know how to set this up without sending the trainPop, parents, chr,
# and useQTL parameters without having them in the Global Environment
# Make a TP nDiffTrainPop times and get the average correlation from that as a
# measure of how good the weights are for the components for optimizing the TP
nDiffTrainPop <- 6

# Scoping problem with parallelization.  The foreach is within the function
# and so I think it creates environments
doThisInParallel <- FALSE
if (doThisInParallel){
maximizeThis <- function(diffLenWgt, ksWgt, covDistWgt){
  require(doParallel); require(foreach)
  cl <- makeCluster(6)
  registerDoParallel(cl)

  allCor <- foreach(i=1:nDiffTrainPop, .combine=c) %dopar% {
    trainOpt <- makeTheTP(diffLenWgt=diffLenWgt,
                          ksWgt=ksWgt,
                          covDistWgt=covDistWgt,
                          candidates=trainCand, parents=parents,
                          tpSize=tpSize, useQTL=useQTL, chr=chr,
                          nIter=1000, verbose=F)
    trainOpt <- AlphaSimR::setPheno(trainOpt, varE=1)
    gblupFromOpt <- genPred(trainOpt, progenyPop)
    cor(gblupFromOpt[tpSize + 1:nProgeny], AlphaSimR::bv(progenyPop))
  }

  stopCluster(cl)
  return(list(Score=mean(allCor, na.rm=T), Pred=0))
}

maximizeThis <- function(diffLenWgt, ksWgt, covDistWgt){
  require(parallel)


  allCor <- foreach(i=1:nDiffTrainPop, .combine=c) %dopar% {
    trainOpt <- makeTheTP(diffLenWgt=diffLenWgt,
                          ksWgt=ksWgt,
                          covDistWgt=covDistWgt,
                          candidates=trainCand, parents=parents,
                          tpSize=tpSize, useQTL=useQTL, chr=chr,
                          nIter=1000, verbose=F)
    trainOpt <- AlphaSimR::setPheno(trainOpt, varE=1)
    gblupFromOpt <- genPred(trainOpt, progenyPop)
    cor(gblupFromOpt[tpSize + 1:nProgeny], AlphaSimR::bv(progenyPop))
  }

  stopCluster(cl)
  return(list(Score=mean(allCor, na.rm=T), Pred=0))
}

} else{

maximizeThis <- function(diffLenWgt, ksWgt, covDistWgt){
  allCor <- NULL
  for (i in 1:nDiffTrainPop){
    trainOpt <- makeTheTP(diffLenWgt=diffLenWgt,
                          ksWgt=ksWgt,
                          covDistWgt=covDistWgt,
                          candidates=trainCand, parents=parents,
                          tpSize=tpSize, useQTL=useQTL, chr=chr,
                          nIter=1000, verbose=F)
    trainOpt <- AlphaSimR::setPheno(trainOpt, varE=1)
    gblupFromOpt <- genPred(trainOpt, progenyPop)
    accFromOpt <- cor(gblupFromOpt[tpSize + 1:nProgeny],
                      AlphaSimR::bv(progenyPop))
    allCor <- c(allCor, accFromOpt)
  }
  return(list(Score=mean(allCor, na.rm=T), Pred=0))
}
}
