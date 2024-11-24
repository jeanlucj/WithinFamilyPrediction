library(rrBLUP)
phenos <- AlphaSimR::setPheno(founders, varE=1, onlyPheno=T)
snpDat <- AlphaSimR::pullSnpGeno(founders) - 1
tstSnpBLUP <- mixed.solve(y=phenos, Z=snpDat, method="REML")

foundersGRM <- A.mat(snpDat)
tstGBLUP <- mixed.solve(y=phenos, K=foundersGRM, method="REML")

snpBLUPfromGBLUP <- t(snpDat) %*% solve(foundersGRM) %*% tstGBLUP$u
cor(tstSnpBLUP$u, snpBLUPfromGBLUP)
