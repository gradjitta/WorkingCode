library("Matrix")
library("gtools")
source("getNbGridsF.R")
source("resetCountsKc.R")
 source("updateCountsS.R")
 source("getTransProbsfcp.R")
source("predHMMdtc.R")
source("resultsScriptn.R")
paths <- readRDS("../grid8train/data/paths8Data.rds")
gridlist <- paths[[4]]
lookup1 <- paths[[5]]
paths1 <- paths[[2]]
cSize <- 8
Lmin <- 100
Bmin <- 80
M <- (floor(Bmin / cSize) +1) * (floor(Lmin / cSize) + 1)
rtt <- readRDS("../grid8train/data/results/resultsalpha3g100C8Train.rds")
tcountsi <- resetCountsKc(M, rtt$K, lookup1, cSize, Lmin, Bmin)
utcounts <- updateCountsS(M, gridlist, paths1, lookup1, rtt$K, rtt$L, tcountsi, cSize)
ptrans <- getTransProbsfcp(M, utcounts, rtt$lw, 10, rtt$K, cSize)
paths1 <- paths[[3]]
startcount <- sapply(1:M, function(i) sum(sapply(1:length(paths1), function(j) paths1[[j]]$grids[1] == i)))
source("resultsScriptnTest.R")
source("predHMMTest.R")
resultslist <- list()
for(i in 1:5) {
  resultslist[[i]] <- resultsScriptnTest(M,100,rtt,i,1,1,startcount,ptrans,cSize,Lmin,Bmin)
}
