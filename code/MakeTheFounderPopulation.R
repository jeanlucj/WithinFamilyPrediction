## ---- Make founder population ------------------------------------------------
# Create haplotypes for founder population of outbred individuals
# Note: default effective population size for runMacs is 100
founderHaps <- AlphaSimR::runMacs(species="WHEAT", nInd=nFounders, nChr=nChr,
                                  segSites=nSnpPerChr+nQTLperChr)

# New global simulation parameters from founder haplotypes
SP <- AlphaSimR::SimParam$new(founderHaps)
SP$restrSegSites(nQTLperChr, nSnpPerChr, overlap=FALSE)
SP$addSnpChip(nSnpPerChr)

# Additive trait architecture
# By default, the genetic variance will be 1
SP$addTraitA(nQtlPerChr=nQTLperChr)

# Create a new population of founders
founders <- AlphaSimR::newPop(founderHaps, simParam=SP)
# Make the founders inbred
founders <- AlphaSimR::makeDH(founders)

## ---- Save important objects -------------------------------------------------
globalEnvFileName <- paste0("saveGlobalEnvObj", random_seed, ".rds")
print(globalEnvFileName)
globalEnvObj <- c("founders", "nChr", "nQTLperChr", "nSnpPerChr",
                  "nFounders", "nTestPop", "nTunePop", "nTrainPop", "nProgeny",
                  "SP")
globalEnvObj <- mget(globalEnvObj) # Get the named objects
saveRDS(globalEnvObj, here::here("output", globalEnvFileName))
