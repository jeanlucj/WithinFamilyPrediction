# Create the founder population that will be the basis for within-family
# prediction The assumption is that different founder populations might behave
# somewhat differently

## ----simulation parameters----------------------------------------------------
# Make sure you set a different random seed each time you want to do this...
random_seed <- 22345
set.seed(random_seed)

nChr <- 7 # Number of chromosomes
nSnpPerChr <- 500 # Number of segregating sites _per chromosome_
nQTLperChr <- 20 # Vary this parameter to get oligo- versus poly- genic traits

nFounders <- 1000
nProgeny <- 100

## ----load packages------------------------------------------------------------
packages_used <- c("AlphaSimR", "tidyverse", "workflowr", "here", "rrBLUP",
                   "EMMREML")
ip <- installed.packages()
all_packages_installed <- TRUE
for (package in packages_used){
  if (!(package %in% ip[,"Package"])){
    print(paste("Please install package", package))
    all_packages_installed <- FALSE
  }
}#END packages_used
if (!all_packages_installed) stop("Need to install packages")

library(tidyverse)

## ----set here path------------------------------------------------------------
here::i_am("code/MakeTheFounderPopulation.R")

## ----Make founder population--------------------------------------------------
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

## ----Save important objects---------------------------------------------------
globalEnvFileName <- paste0("saveGlobalEnvObj", random_seed, ".rds")
globalEnvObj <- c("founders", "nChr", "nFounders", "nProgeny", "nQTLperChr",
                  "nSnpPerChr", "SP")
globalEnvObj <- mget(globalEnvObj)
saveRDS(globalEnvObj, here::here("output", globalEnvFileName))
