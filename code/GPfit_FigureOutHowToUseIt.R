n = 5
d = 1
computer_simulator <- function(x){
  x = 2 * x + 0.5
  y = sin(10 * pi * x) / (2 * x) + (x - 1)^4
  return(y)
}
set.seed(3)
library(lhs)
x = maximinLHS(n, d)
y = computer_simulator(x)
GPmodel = GP_fit(x, y, corr=list(type = "matern", nu=5/2))
GPmodel = GP_fit(x, y, nug_thres=10)
GPmodel = GP_fit(x, y)
print(GPmodel)

pltpts <- seq(from=0, to=1, length.out=11)
plot(pltpts, computer_simulator(pltpts), pch=16)
points(pltpts, predict(GPmodel, xnew=pltpts)$Y_hat, pch=16, col="blue")

stdUnit <- function(v) (v - min(v)) / (range(v) %>% diff)
initObsStd <- apply(initializeObs[,-c(1, 5)], 2, stdUnit)
gpIO <- GP_fit(initObsStd[1:50,], unlist(initializeObs[1:50, 5]),
               corr=list(type="matern", nu=5/2),
               nug_thres = 2)
print(gpIO)

nObs <- 100
gpRes <- NULL
for (nt in c(0.5, 1, 1.5, 2, 2.5, 3, 10)){
  gpIO <- GP_fit(initObsStd[1:nObs,], unlist(initializeObs[1:nObs, 5]),
                 corr=list(type="matern", nu=5/2),
                 nug_thres = nt)
  gpRes <- rbind(gpRes, c(nt, gpIO$beta, gpIO$sig2, gpIO$delta))
}


zeroToOne <- seq(from=0, to=1, length.out=11)
gridHyperCube <- cbind(rep(zeroToOne, each=121),
                       rep(rep(zeroToOne, each=11), times=11),
                       rep(zeroToOne, times=121))

nt <- 8
nObs <- nrow(initializeObs)
gpIO <- GP_fit(initObsStd[1:nObs,], unlist(initializeObs[1:nObs, 5]),
               corr=list(type="matern", nu=5/2),
               nug_thres = nt)
print(gpIO)
tstMat <- predict(gpIO, xnew=gridHyperCube)
tstMat$complete_data[which.max(tstMat$Y_hat),]
