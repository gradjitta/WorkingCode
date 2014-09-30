require(snowfall)
#sfSetMaxCPUs(number=80)
paths <- readRDS("Data/paths.rds")
L <- readRDS("newData/inLc.rds")
lookup <- readRDS("newData/lookup1.rds")
new1000 <- readRDS("newData/new1000c.rds")
tcounts <- readRDS("newData/tcountsc.rds")
states <- new1000$CK
M <- 520
N <- length(paths)
bookL <- readRDS("newData/bookLc.rds")
source("SampleDirichlet.R")
source("generateV.R")
source("alphaInit.R")
source("getTransp.R")
source("getTransProbs.R")
source("getNbGrids.R")
source("getEmitf.R")
source("pFFBSIhmmf.R")
source("updateHDP.R")
source("updateLpi.R")
source("sampleVk.R")
source("getWeights.R")
source("sampleAjk1.R")
source("sampleMeanjnm.R")
source("sampleVariancejn.R")
source("lprobGImp.R")
source("isampleAjk1.R")
source("resetCounts1.R")
source("dataForHdp.R")
library(Matrix)
iterNum <- 3
#sfInit(socketHosts=rep(c('ukko156', 'ukko099', 'ukko153', 'ukko101', 'ukko102', 'ukko109', 'ukko024', 'ukko100', 'ukko143', 'ukko144'),each=8), cpus=80,type='SOCK',parallel=T)
#sfInit(socketHosts=rep(c('ukko156','ukko099','ukko153','ukko101','ukko102','ukko109','ukko100','ukko024'),each=8), cpus=64,type='SOCK',parallel=T)
#sfInit(socketHosts=rep(c('ukko102','ukko109','ukko100','ukko024','ukko023'),each=8), cpus=40,type='SOCK',parallel=T)
sfInit(socketHosts=rep(c('ukko102','ukko109','ukko100'),each=8), cpus=24,type='SOCK',parallel=T)
sfExport("paths","lookup")
print("paths exported..")
sfExport("SampleDirichlet","alphaInit","getTransp","getNbGrids","getEmitf","pFFBSIhmmf","updateHDP",
"updateLpi","sampleVk","getWeights","sampleAjk1","sampleMeanjnm","sampleVariancejn","lprobGImp","isampleAjk1","resetCounts1","dataForHdp", "getTransProbs")
print("functions exported..")
N <- length(paths)
iter <- 0
while(iter < iterNum) {
if(iter > 0) {
   print("Reset counts..")
   tcounts <- resetCounts1(M, states, lookup)
   #sfExport("paths","tcounts","L","new1000","lookup")
   #print(new1000)
}
#sfInit(socketHosts=rep(c('ukko102','ukko109','ukko100','ukko024','ukko023'),each=8), cpus=40,type='SOCK',parallel=T)
print(sum(all.equal(tcounts, resetCounts1(M, states, lookup)) == TRUE))
#tcounts <- resetCounts1(M, states, lookup)
#sfExport("pFFBSIhmm","getTrans", "getEmit", "getNbGrids","SampleDirichlet", "alphaInit","sampleVk","getWeights","")
#sfExport("SampleDirichlet.R","alphaInit.R","getTrans.R","getNbGrids.R","getEmit.R","pFFBSIhmm.R","updateHDP.R","updateLpi.R","sampleVk.R","getWeights.R","sampleAjk1.R","sampleMeanjnm.R","sampleVariancejn.R","lprobGImp.R","isampleAjk1.R","resetCounts1.R","dataForHdp.R")
time1 <- Sys.time()
tprobs <- getTransProbs(tcounts)
sfExport("tprobs","L","new1000")
print("before pffbs..")
temps.Flapply <- sfLapply(1:length(paths), pFFBSIhmmf, paths, tprobs, L, new1000, lookup)
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
sfStop()
