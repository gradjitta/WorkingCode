set.seed(4673)
iterNum <- 50
paths <- readRDS("newData/paths.rds")
L <- readRDS("newData/inLc.rds")
lookup <- readRDS("newData/lookup1.rds")
new1000 <- readRDS("newData/new1000c.rds")
tcounts <- readRDS("newData/tcountsc.rds")
spaths <- readRDS("newData/spaths.rds")
spaths1 <- spaths[[1]]
spaths2 <- spaths[[2]]
source("SampleDirichlet.R")
source("generateV.R")
source("alphaInit.R")
source("getTransp.R")
source("getNbGrids.R")
source("getEmitf.R")
source("pFFBSIhmmf.R")
source("updateHDP.R")
source("updateLpic.R")
source("sampleVk.R")
source("getWeights.R")
source("sampleAjk1.R")
source("sampleMeanjnm.R")
source("sampleVariancejn.R")
source("lprobGImp.R")
source("isampleAjk1.R")
source("resetCounts1.R")
source("dataForHdp.R")
source("collectCounts.R")
source("updateCounts2.R")
source("getTransProbs.R")
library(Matrix)
require(snowfall)
library(parallel)
# Initialize transition data structure
states <- new1000$CK
d <- new1000$d
phi <- new1000$phi
bk <- new1000$beta.k
G <- length(bk)
# Main Loop starts
sfInit(socketHosts=rep(c('ukko102','ukko109','ukko100'),each=8), cpus=24,type='SOCK',parallel=T)
sfExport("paths","lookup")
print("paths exported..")
sfExport("SampleDirichlet","alphaInit","getTransp","getNbGrids","getEmitf","pFFBSIhmmf","updateHDP",
"updateLpic","sampleVk","getWeights","sampleAjk1","sampleMeanjnm","sampleVariancejn","lprobGImp","isampleAjk1","resetCounts1","dataForHdp", "getTransProbs")
print("functions exported..")
N <- length(paths)
M <- 520
########################3
iter <- 0
while(iter < iterNum) {
   if(iter > 0) {
      print("Reset counts..")
      tcounts <- resetCounts1(M, states, lookup)
   }
   time4 <- Sys.time()
   #print(sum(all.equal(tcounts, resetCounts1(M, states, lookup)) == TRUE))
   print(iter)
   print("update counts")
   tempCounts <- mclapply(1:10, updateCounts2, 520, paths, lookup, new1000$CK, tcounts, L, spaths1, spaths2, mc.cores = 10)
   print("combine counts")
   tcountsAgg <- collectCounts(tempCounts)
   tcounts <- list(tcountsAgg[[1]],tcountsAgg[[2]],tcountsAgg[[3]])
   bookL <- tcountsAgg[[4]]
   print(Sys.time()- time4)
   time1 <- Sys.time()
   tprobs <- getTransProbs(tcounts)
   sfExport("tprobs","L","new1000")
   print("before pffbs..")
   L <- sfLapply(1:length(paths), pFFBSIhmmf, paths, tprobs, L, new1000, lookup)
   sfRemove("tprobs","L","new1000")
   print("pFFBS done..")
   print(Sys.time() - time1)
   # 6) Sample phi, local sticks similar to HDP
   # L (but d for sample phi) from backward sampling used for the HDP updates or D
   Doom <- dataForHdp(paths, L, states, lookup)
   data <- Doom$data
   Lhdp <- Doom$Lhdp
   d <- Doom$d
   updates <- updateHDP(data, states, Lhdp, d, lookup, tcounts, bookL)
   bk <- updates$beta.k
   states <- updates$CK
   phi <- updates$phi
   lpi <- updates$lpi
   new1000$CK <- states
   new1000$phi <- phi
   new1000$d <- d
   new1000$beta.k <- bk
   new1000$lpi <- lpi
   new1000$L <- Lhdp
   iter <- iter + 1
}
saveRDS(tcounts,"results/tcounts50.rds")
saveRDS(L,"results/Lf50.rds")
saveRDS(new1000, "results/new10050.rds")
sfStop()
