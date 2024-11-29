# Create the founder population that will be the basis for within-family
# prediction The assumption is that different founder populations might behave
# somewhat differently

## ----simulation parameters----------------------------------------------------
nChr <- 7 # Number of chromosomes
nSnpPerChr <- 500 # Number of markers _per chromosome_
nQTLperChr <- 20 # Vary this parameter to get oligo- versus poly- genic traits

nFounders <- 1000
nTestPop <- 200
nTunePop <- 200
nTrainPop <- nFounders - nTestPop - nTunePop
nProgeny <- 100 # Number of progeny in test biparental

## ----load packages------------------------------------------------------------
packages_used <- c("AlphaSimR",
                   "tidyverse", "tidymodels",
                   "workflowr", "here",
                   "rrBLUP", "EMMREML")
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
here::i_am("code/SetSimulationParameters.R")

print(utils::sessionInfo())
