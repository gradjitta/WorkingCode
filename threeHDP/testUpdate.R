require(snowfall)
library(Matrix)
#sfSetMaxCPUs(number=80)
paths <- readRDS("Data/paths.rds")
L <- readRDS("Data/inL.rds")
lookup <- readRDS("newData/lookup1.rds")
new1000 <- readRDS("newData/new1000.rds")
tcounts <- readRDS("newData/tcountsnew.rds")
Doom <- readRDS("Data/doom.rds")
states <- Doom$Ck
res <- list()
M <- 520
N <- length(paths)
source("resetCounts1.R")
source("updateCounts2.R")
source("getNbGrids.R")
initcounts <- resetCounts1(M, states, lookup)
sfInit(socketHosts=rep(c('ukko099','ukko044','ukko098','ukko175','ukko144'),each=8), cpus=40,type='SOCK',parallel=T)
sfLibrary(Matrix)
sfExport("M","paths","lookup","states","initcounts","L")
sfExport("getNbGrids","resetCounts1", "updateCounts2")
res <- sfLapply(1:length(paths), updateCounts2, M, paths, lookup, states, initcounts, L)
t1 <- res[[1]][[1]]
t2 <- res[[1]][[2]]
tg <- res[[1]][[3]]
for (p in 2:N) {
   for (i in 1:520) {
      J <- length(res[[p]][[1]][[i]])
      for (j in 1:J) {
         if (!is.null(res[[p]][[1]][[i]][[j]])) {
           t1[[i]][[j]] <- t1[[i]][[j]] + res[[p]][[1]][[i]][[j]]
         }  
      }
   }
}
for (p in 2:N) {
  tg <- tg + res[[p]][[3]]
}
saveRDS(t1,"results/tempt1.rds")
saveRDS(tg,"results/tg.rds")
sfStop()
