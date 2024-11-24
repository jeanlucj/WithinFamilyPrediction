# Get a baseline accuracy for a random training population
here::i_am("code/AccuracyOfRandomTPs.R")
globalEnvObj <- readRDS(here::here("output", "saveGlobalEnvObj.rds"))
for (n in names(globalEnvObj)) assign(n, globalEnvObj[[n]], envir=.GlobalEnv)
allRes <- NULL
for (tpSize in c(50, 100, 200, 400)){
  tpSzRes <- NULL
  for (i in 1:400){
    chr <- 1

    require(tidyverse)
    source(here::here("code", "functionsToUseBO.R"))
    source(here::here("code", "genomicPredictionFunction.R"))
    source(here::here("code", "assessIndMrkSimilarityToPar.R"))
    source(here::here("code", "markerCovarianceDistance.R"))

    # Make the progeny pop
    parentsIdx <- sample(nFounders, 2)
    parents <- founders[parentsIdx]
    trainCand <- founders[-parentsIdx]
    f1 <- AlphaSimR::makeCross(parents, crossPlan=matrix(1:2, nrow=1))
    progenyPop <- AlphaSimR::makeDH(f1, nDH=nProgeny)

    trainRnd <- trainCand[sample(AlphaSimR::nInd(trainCand), tpSize)]
    trainRnd <- AlphaSimR::setPheno(trainRnd, varE=1)
    gblupFromRnd <- genPred(trainRnd, progenyPop)[tpSize + 1:nProgeny]
    tpSzRes <- c(tpSzRes,
                 cor(gblupFromRnd, AlphaSimR::bv(progenyPop)))
    if (i %% 20 == 0) cat(".")
  }
  cat("\n")
  allRes <- cbind(allRes, tpSzRes)
}
colnames(allRes) <- paste0("tpSz", c(050, 100, 200, 400))
saveRDS(allRes, here::here("output", "AccuracyOfRandomTPs.rds"))
apply(allRes, 2, summary, na.rm=T)
boxplot(allRes)
