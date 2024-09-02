# Use RRBLUP to run genomic prediction
genPred <- function(traiPop, testPop){
  require(rrBLUP)
  nTest <- AlphaSimR::nInd(testPop)
  df <- data.frame(phenoVal=c(traiPop@pheno, rep(NA, nTest)),
                   genoID=c(traiPop@id, testPop@id))
  snp <- rbind(AlphaSimR::pullSnpGeno(traiPop),
               AlphaSimR::pullSnpGeno(testPop)) - 1
  grm <- rrBLUP::A.mat(snp)
  gBLUP <- rrBLUP::kin.blup(df, geno="genoID", pheno="phenoVal", K=grm)
  return(gBLUP$g)
}
