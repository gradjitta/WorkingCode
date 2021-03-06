paths <- readRDS("newData/paths.rds")
L <- readRDS("newData/inLc.rds")
lookup <- readRDS("newData/lookup1.rds")
tcountsi <- readRDS("newData/tcountsINIT.rds")
gridlist <- readRDS("newData/gridlist.rds")
source("getTransProbsf.R")
source("updateCountsNew.R")
source("pFFBSIhmmf.R")
source("alphaInit.R")
source("getTransp.R")
source("getEmitf.R")
source("getNbGrids.R")
library(Matrix)
require(MCMCpack)
require(snowfall)
library(parallel)
new1000 <- readRDS("newData/new1000c.rds")
print("tcounts -..")
tcounts <- updateCountsNew(520, gridlist, paths, lookup, new1000$CK, L, tcountsi)
tprobs <- getTransProbsf(tcounts)
print("main algo..")
time1 <- Sys.time()
mclapply(1:length(paths), pFFBSIhmmf, paths,tprobs, L, new1000, lookup, mc.cores = 12)
print(Sys.time() - time1)
