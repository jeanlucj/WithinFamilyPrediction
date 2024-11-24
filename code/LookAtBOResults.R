here::i_am("code/LookAtBOResults.R")

initializeObs <- readRDS(here::here("output", "ParmsOptimizedByBOparallel.rds"))
initializeObs <- initializeObs$History

plot(initializeObs$diffLenWgt, pch=16, ylab="diffLenWgt")
plot(initializeObs$ksWgt, pch=16, ylab="ksWgt")
plot(initializeObs$covDistWgt, pch=16, ylab="covDistWgt")
plot(initializeObs$Value, pch=16, ylab="Accuracy")

fitCDW <- lm(covDistWgt ~ Round, data=initializeObs)
anova(fitCDW)
summary(fitCDW)

fitAcc <- lm(Value ~ Round, data=initializeObs)
anova(fitAcc)
summary(fitAcc)

# Compare the variance in the first nStrt versus the last nEnd of the obs
nS <- 50
nE <- 50
strtEndSD <- function(v, nStrt=nS, nEnd=nE){
  lenV <- length(v)
  return(round(c(sd(v[1:nStrt], na.rm=T), sd(v[1:nEnd + lenV - nEnd], na.rm=T)), 3))
}
apply(initializeObs[,-1], 2, strtEndSD)
boxplot(initializeObs$Value[(1 - nE):0 + nrow(initializeObs)])
apply(initializeObs[(1 - nE):0 + nrow(initializeObs),-1], 2, summary)
apply(10^initializeObs[(1 - nE):0 + nrow(initializeObs), c(-1, -5)], 2, summary) %>% round(2)
# cor(initializeObs[(1 - nE):0 + nrow(initializeObs),-1])

# Trying to use some local smoothing with GPfit or Loess to find a maximum
# It's not giving results that I believe
apply(initializeObs[,-1], 2, summary)
obsHyperCube <- apply(initializeObs[,2:4], 2, function(v) 
  return((v - min(v)) / (range(v) %>% diff)))
GPmodelMat <- GPfit::GP_fit(obsHyperCube, initializeObs[,5], 
                            corr=list(type = "matern", nu=5/2))
zeroToOne <- seq(from=0, to=1, length.out=11)
gridHyperCube <- cbind(rep(zeroToOne, each=121),
                       rep(rep(zeroToOne, each=11), times=11),
                       rep(zeroToOne, times=121))
tstMat <- predict(GPmodelMat, xnew=gridHyperCube)
tstMat$complete_data[which.max(tstMat$Y_hat),]

gridHyperCube2 <- data.frame((gridHyperCube - 0.5) * 3.4)
colnames(gridHyperCube2) <- c("diffLenWgt", "ksWgt", "covDistWgt")
loessFit <- loess(Value ~ diffLenWgt + ksWgt + covDistWgt, 
                  data=initializeObs)
tstLoess <- predict(loessFit, newdata=gridHyperCube2)
gridHyperCube2[which.max(tstLoess),]

tstLoessObsDat <- predict(loessFit)
c(unlist(initializeObs[which.max(tstLoessObsDat),]), max(tstLoessObsDat))
