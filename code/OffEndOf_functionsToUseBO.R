
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

} else{ # Don't run TP optimization in parallel

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
