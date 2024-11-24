tst <- readRDS(here::here("output", "IndependenceOfOptimizedTPs_FourCrit100.rds"))
t3 <- colMeans(tst, na.rm=T)
xRange <- range(c(100, t3[0:14*4+1], t3[0:14*4+86])) # Still working on this
yRange <- range(c(t3[c(61, 146)], t3[0:14*4+2], t3[0:14*4+4], t3[0:14*4+87], t3[0:14*4+89]))
plot(c(100, t3[0:14*4+1]), c(t3[61], t3[0:14*4+2]),
     pch=16, xlim=xRange, ylim=yRange,
     xlab="Mean N Training Individuals",
     ylab="Mean Accuracy",
     main="FourCrit100")
points(c(100, t3[0:14*4+1]), c(t3[61], t3[0:14*4+4]), pch=16, col="dark green")
points(c(100, t3[0:14*4+86]), c(t3[146], t3[0:14*4+87]), pch=16, col="blue")
points(c(100, t3[0:14*4+86]), c(t3[146], t3[0:14*4+89]), pch=16, col="dark red")

tst <- readRDS(here::here("output", "IndependenceOfOptimizedTPs_FourCrit50.rds"))
t3 <- colMeans(tst)
xRange <- range(c(50, t3[0:14*3+1], t3[0:14*3+71]))
yRange <- range(c(t3[c(46, 116)], t3[0:14*3+2], t3[0:14*3+72]))
plot(c(50, t3[0:14*3+1]), c(t3[46], t3[0:14*3+2]),
     pch=16, xlim=xRange, ylim=yRange,
     xlab="Mean N Training Individuals",
     ylab="Mean Accuracy")
points(c(50, t3[0:14*3+71]), c(t3[116], t3[0:14*3+72]), pch=16, col="blue")

tst <- readRDS(here::here("output", "IndependenceOfOptimizedTPs_FourCrit.rds"))
t3 <- colMeans(tst)
xRange <- range(c(100, t3[0:6*3+1], t3[0:6*3+35]))
yRange <- range(c(t3[c(22, 56)], t3[0:6*3+2], t3[0:6*3+36]))
plot(c(100, t3[0:6*3+1]), c(t3[22], t3[0:6*3+2]),
     pch=16, xlim=xRange, ylim=yRange,
     xlab="Mean N Training Individuals",
     ylab="Mean Accuracy")
points(c(100, t3[0:6*3+35]), c(t3[56], t3[0:6*3+36]), pch=16, col="blue")

t3median <- apply(tst, 2, median)
xRange <- range(c(100, t3median[0:6*3+1], t3median[0:6*3+35]))
yRange <- range(c(t3median[c(22, 56)], t3median[0:6*3+2], t3median[0:6*3+36]))
plot(c(100, t3median[0:6*3+1]),
     c(t3median[22], t3median[0:6*3+2]),
     pch=16, xlim=xRange, ylim=yRange,
     xlab="Median N Training Individuals",
     ylab="Median Accuracy")
points(c(100, t3median[0:6*3+35]), c(t3median[56], t3median[0:6*3+36]), pch=16, col="blue")

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
