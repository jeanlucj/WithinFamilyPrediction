# Use RRBLUP to run genomic prediction
genPred <- function(trainPop, testPop){
  require(rrBLUP)
  nTest <- AlphaSimR::nInd(testPop)
  df <- data.frame(phenoVal=c(trainPop@pheno, rep(NA, nTest)),
                   genoID=c(trainPop@id, testPop@id))
  snp <- rbind(AlphaSimR::pullSnpGeno(trainPop),
               AlphaSimR::pullSnpGeno(testPop)) - 1
  grm <- rrBLUP::A.mat(snp)
  gBLUP <- rrBLUP::kin.blup(df, geno="genoID", pheno="phenoVal", K=grm)
  return(gBLUP$g)
}

# For the testPop, you only want to predict the BLUPs for the chromosome
# allTrain population with all individuals that have phenotypes
# chrTrain index of the individuals in that population to train for the chromo
# chr what chromosome are you training for
# The version calculates ONE kernel / variance component for the focal chromo
# It would be good to have a version with one kernel for chrTrain and another
# for the remainder of all training individuals. That way you could have two
# variance components. (Why is that a good thing?).  It's not that I want two
# variance components, but I want the model to avoid having to assume that there
# is zero covariance between chrTrain and restTrain.  Doesn't that get assumed
# anyway? Check out the emmremlMultiKernel documentation. Yes, it does.
genPredOneChr <- function(allTrain, chrTrain, testPop, chr){
  require(EMMREML); require(rrBLUP)
  nTrain <- AlphaSimR::nInd(allTrain)
  nTest <- AlphaSimR::nInd(testPop)
  # Split the SNPs into focal Chromosome (fChr) versus rest of genome (Rog)
  snp <- rbind(AlphaSimR::pullSnpGeno(allTrain),
               AlphaSimR::pullSnpGeno(testPop)) - 1
  snpNames <- colnames(snp)
  sn <- strsplit(snpNames, "_")
  fChr <- which(sapply(sn, function(v) v[1]) == as.character(chr))
  Rog <- setdiff(1:ncol(snp), fChr)
  # Make GRMs with the fChr and Rog snps
  grmChr <- rrBLUP::A.mat(snp[, fChr])
  grmRog <- rrBLUP::A.mat(snp[, Rog])
  # You want the chrTrain individuals to share information with the testPop for
  # the focal chromosome but _not_ the other individuals in allTrain
  # That means setting to zero all covariances that are between the groups
  # {chrTrain, testPop} and setdiff(alltrain, {chrTrain, testPop}) in the grmChr
  chrTandTest <- c(chrTrain, nTrain + 1:nTest)
  nonChrT <- setdiff(1:(nTrain+nTest), chrTandTest)
  grmChr[chrTandTest, nonChrT] <- 0
  grmChr[nonChrT, chrTandTest] <- 0
  phenos <- allTrain@pheno
  fixedEff <- matrix(rep(1, nTrain), ncol=1)
  zChr <- zRog <- diag(nrow=nTrain, ncol=nTrain+nTest)
  mmRes <- EMMREML::emmremlMultiKernel(y=phenos, X=fixedEff,
                              Zlist=list(zChr, zRog),
                              Klist=list(grmChr, grmRog))
  return(mmRes$uhat)
}

# Here, it's assumed that you enter with a list of training populations that has
# as many training populations as there are chromosomes
genPredMultiChr <- function(allTrain, trainPopList, testPop){
  require(EMMREML); require(rrBLUP)

  nTrain <- AlphaSimR::nInd(allTrain)
  nTest <- AlphaSimR::nInd(testPop)
  allChrPred <- numeric(nTest)
  nChr <- length(trainPopList)
  for (chr in 1:nChr){
    chrPred <- genPredOneChr(allTrain, trainPopList[[chr]], testPop, chr)
    allChrPred <- allChrPred + chrPred[nTrain + 1:nTest]
  }
  return(allChrPred)
}

# There are two SNP matrices, one for portions of the genome that are presumed
# IBD with the parents and one for portions that are presumed not IBD
# This function returns the SNP effects from the ibd matrix
genPredPartitionedSnps <- function(ibdMat, nonIbdMat, phenos){
  require(EMMREML); require(rrBLUP)
  ibdGRM <- rrBLUP::A.mat(ibdMat*2 - 1)
  nonIbdGRM <- rrBLUP::A.mat(nonIbdMat*2 - 1)
  nInd <- length(phenos)
  fixedEff <- matrix(rep(1, nInd), ncol=1)
  zIbd <- zNonIbd <- diag(nrow=nInd)
  mmRes <- EMMREML::emmremlMultiKernel(y=phenos, X=fixedEff,
                                       Zlist=list(zIbd, zNonIbd),
                                       Klist=list(ibdGRM, nonIbdGRM))
  ibdSnpEff <- t(ibdMat) %*% solve(ibdGRM + diag(1e-5, nInd)) %*% mmRes$uhat[1:nInd]
  return(ibdSnpEff)
}
