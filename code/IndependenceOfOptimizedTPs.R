# Independence of separately optimized TPs
here::i_am("code/IndependenceOfOptimizedTPs.R")
# Make progeny populations, optimize the TP five times, find out how many
# candidates are covered
allRes <- NULL
for (parPair in 1:100){ # Leave this in the loop in case I parallelize it
  globalEnvObj <- readRDS(here::here("output", "saveGlobalEnvObj.rds"))
  for (n in names(globalEnvObj)) assign(n, globalEnvObj[[n]], envir=.GlobalEnv)
  # For now using the criterion weights that came from ~280 iterations of BO
  # These values used by mistake: only covDist mattered then because 10^wgt
  # diffLenWgt=0.97; ksWgt=0.32; covDistWgt=24.19
  diffLenWgt=-0.5; ksWgt=1; covDistWgt=1.3
  tpSize <- 100
  useQTL <- F
  chr <- 1
  nIterTpOptim <- 1000
  verbose <- T
  chr <- 1

  # How many training populations to test per family
  nTPperFam <- 24

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

  # First do multiple TPs for optimized TPs
  perTPAcc <- NULL
  collectTPids <- NULL
  collectPreds <- NULL
  collectStdPreds <- NULL

  for (i in 1:nTPperFam){
    trainOpt <- makeTheTP(diffLenWgt=diffLenWgt,
                          ksWgt=ksWgt,
                          covDistWgt=covDistWgt,
                          candidates=trainCand, parents=parents,
                          tpSize=tpSize, useQTL=useQTL, chr=chr,
                          nIter=nIterTpOptim, verbose=verbose)
    trainOpt <- AlphaSimR::setPheno(trainOpt, varE=1)
    gblupFromOpt <- genPred(trainOpt, progenyPop)[tpSize + 1:nProgeny]
    perTPAcc <- c(perTPAcc,
                  cor(gblupFromOpt, AlphaSimR::bv(progenyPop)))
    if (i == 1){
      collectTPids <- trainOpt@iid
    } else{
      collectTPids <- cbind(collectTPids, trainOpt@iid)
    }
    collectPreds <- cbind(collectPreds, gblupFromOpt)
    collectStdPreds <- cbind(collectStdPreds, scale(gblupFromOpt))
  }# END optimize separate TPs
  meanAccPerTP <- mean(perTPAcc)

  # Maybe it's worth making a lot of independently optimized TPs
  # How fast does it saturate?
  nSeparateTPs <- 16

  collectNTPAcc <- NULL
  for (nTPs in 2:nSeparateTPs){
    nRepSeparateTPs <- 30
    allNTP <- allAccMeanPreds <- allAccMeanStdPreds <- allUniTPAcc <- NULL
    for (i in 1:nRepSeparateTPs){
      tpIdx <- sample(nTPperFam, nTPs)
      trainUni <- collectTPids[,tpIdx] %>% c %>% unique
      nTPids <- trainUni %>% length
      accMeanPreds <- cor(apply(collectPreds[,tpIdx], 1, mean),
                          AlphaSimR::bv(progenyPop))
      accMeanStdPreds <- cor(apply(collectStdPreds[,tpIdx], 1, mean),
                             AlphaSimR::bv(progenyPop))
      allNTP <- c(allNTP, nTPids)
      allAccMeanPreds <- c(allAccMeanPreds, accMeanPreds)
      allAccMeanStdPreds <- c(allAccMeanStdPreds, accMeanStdPreds)
      # What accuracy do you get if you train on these individuals?
      trainUni <- trainCand[trainCand@iid %in% trainUni]
      trainUni <- AlphaSimR::setPheno(trainUni, varE=1)
      gblupFromUni <- genPred(trainUni, progenyPop)[nTPids + 1:nProgeny]
      allUniTPAcc <- c(allUniTPAcc,
                       cor(gblupFromUni, AlphaSimR::bv(progenyPop)))

    }
    collectNTPAcc <- c(collectNTPAcc,
                       mean(allNTP),
                       mean(allAccMeanPreds),
                       mean(allAccMeanStdPreds),
                       mean(allUniTPAcc))
  }#END 1:nSeparateTP
  # allNTPAcc is the big one
  allNTPAcc <- c(collectNTPAcc, meanAccPerTP, perTPAcc)

  # Second do this for random TPs
  perTPAcc <- NULL
  collectTPids <- NULL
  collectPreds <- NULL
  collectStdPreds <- NULL

  for (i in 1:nTPperFam){
    trainRnd <- trainCand[sample(AlphaSimR::nInd(trainCand), tpSize)]
    trainRnd <- AlphaSimR::setPheno(trainRnd, varE=1)
    gblupFromRnd <- genPred(trainRnd, progenyPop)[tpSize + 1:nProgeny]
    perTPAcc <- c(perTPAcc,
                  cor(gblupFromRnd, AlphaSimR::bv(progenyPop)))
    if (i == 1){
      collectTPids <- trainRnd@iid
    } else{
      collectTPids <- cbind(collectTPids, trainRnd@iid)
    }
    collectPreds <- cbind(collectPreds, gblupFromRnd)
    collectStdPreds <- cbind(collectStdPreds, scale(gblupFromRnd))
  }# END random separate TPs
  meanAccPerTP <- mean(perTPAcc)

  collectNTPAcc <- NULL
  for (nTPs in 2:nSeparateTPs){
    nRepSeparateTPs <- 30
    allNTP <- allAccMeanPreds <- allAccMeanStdPreds <- allUniTPAcc <- NULL
    for (i in 1:nRepSeparateTPs){
      tpIdx <- sample(nTPperFam, nTPs)
      trainUni <- collectTPids[,tpIdx] %>% c %>% unique
      nTPids <- trainUni %>% length
      accMeanPreds <- cor(apply(collectPreds[,tpIdx], 1, mean),
                          AlphaSimR::bv(progenyPop))
      accMeanStdPreds <- cor(apply(collectStdPreds[,tpIdx], 1, mean),
                             AlphaSimR::bv(progenyPop))
      allNTP <- c(allNTP, nTPids)
      allAccMeanPreds <- c(allAccMeanPreds, accMeanPreds)
      allAccMeanStdPreds <- c(allAccMeanStdPreds, accMeanStdPreds)
      # What accuracy do you get if you train on these individuals?
      trainUni <- trainCand[trainCand@iid %in% trainUni]
      trainUni <- AlphaSimR::setPheno(trainUni, varE=1)
      gblupFromUni <- genPred(trainUni, progenyPop)[nTPids + 1:nProgeny]
      allUniTPAcc <- c(allUniTPAcc,
                       cor(gblupFromUni, AlphaSimR::bv(progenyPop)))
    }
    collectNTPAcc <- c(collectNTPAcc,
                       mean(allNTP),
                       mean(allAccMeanPreds),
                       mean(allAccMeanStdPreds),
                       mean(allUniTPAcc))
  }#END 1:nSeparateTP
  allNTPAcc <- c(allNTPAcc, collectNTPAcc, meanAccPerTP, perTPAcc)

  allRes <- rbind(allRes, allNTPAcc)

  saveRDS(allRes, here::here("output", "IndependenceOfOptimizedTPs_FourCrit100.rds"))
}#END choose a new parent pair

# Make the column names for allRes
cnAllRes <- NULL
for (nTPs in 2:nSeparateTPs){
  cnAllRes <- c(cnAllRes, paste0(c("nUniqueInd_Opt",
                                   "accMeanPreds_Opt",
                                   "accMeanStdPreds_Opt",
                                   "accMeanUni_Opt"), nTPs))
}
cnAllRes <- c(cnAllRes, "meanAccPerTP_Opt", paste0("perTPAcc_Opt", 1:nTPperFam))

for (nTPs in 2:nSeparateTPs){
  cnAllRes <- c(cnAllRes, paste0(c("nUniqueInd_Rnd",
                                   "accMeanPreds_Rnd",
                                   "accMeanStdPreds_Rnd",
                                   "accMeanUni_Rnd"), nTPs))
}
cnAllRes <- c(cnAllRes, "meanAccPerTP_Rnd", paste0("perTPAcc_Rnd", 1:nTPperFam))

colnames(allRes) <- cnAllRes
saveRDS(allRes, here::here("output", "IndependenceOfOptimizedTPs_FourCrit100.rds"))
