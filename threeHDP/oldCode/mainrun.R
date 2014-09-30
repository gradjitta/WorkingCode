set.seed(4673)
numIter <- 5
paths <- readRDS("Data/paths.rds")
L <- readRDS("newData/inLc.rds")
lookup <- readRDS("newData/lookup1.rds")
new1000 <- readRDS("newData/new1000c.rds")
source("SampleDirichlet.R")
source("generateV.R")
source("alphaInit.R")
source("getTrans.R")
source("getNbGrids.R")
source("getEmit.R")
source("pFFBSIhmm.R")
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
require(snowfall)
#sfSetMaxCPUs(number=80)
#sfInit(socketHosts=rep(c('ukko156', 'ukko101', 'ukko153', 'ukko102', 'ukko109', 'ukko100', 'ukko024', 'ukko168', 'ukko098', 'ukko099'),each=8), cpus=80,type='SOCK',parallel=T)
#sfExport("getTrans", "getEmit", "getNbGrids","SampleDirichlet", "alphaInit")
#sfExportAll()
#print("exports all..")
N <- length(paths)
M <-520
sfInit(socketHosts=rep(c('ukko102','ukko109','ukko100','ukko098','ukko168','ukko044','ukko099'),each=8), cpus=56,type='SOCK',parallel=T)
sfExport("paths","lookup")
print("paths exported..")
sfExport("SampleDirichlet","alphaInit","getTrans","getNbGrids","getEmit","pFFBSIhmm","updateHDP",
"updateLpi","sampleVk","getWeights","sampleAjk1","sampleMeanjnm","sampleVariancejn","lprobGImp","isampleAjk1","resetCounts1","dataForHdp")
print("functions exported..")

# Initialize transition data structure
states <- new1000$CK
d <- new1000$d
phi <- new1000$phi
bk <- new1000$beta.k
G <- length(bk)
tcounts1 <- list() # P(L(t) | L(t-1), S(t))
tcounts2 <- list() # P(S(t) | S(t-1), L(t-1))
tcountsg <- Matrix(0, ncol= M, nrow = M) # for Jumps/ Missing transitions
corners <- c(1,26,495,520)
tcounts1 <- lapply(1:M, function(i) lapply(1:length(getNbGrids(i)), function(j) NULL ))
tcounts2 <- lapply(1:M, function(i) NULL)
for (i in 1:M) {
  if (sum(corners == i) == 0) {
    near.i <- getNbGrids(i)
    for(j in 1:length(near.i)) {
      ti <- which(lookup == i)
      tj <- which(lookup == near.i[j])
      if(sum(ti) != 0 && sum(tj) != 0) {
        tcounts1[[i]][[j]] <- Matrix(0, ncol = length(states[[tj]]),
                                     nrow = length(states[[ti]]))
      } else {
        if(sum(ti) == 0 && sum(tj) != 0) {
          tcounts1[[i]][[j]] <- Matrix(0, ncol = length(states[[tj]]), nrow = 1)
        } else if(sum(ti) != 0 && sum(tj) == 0){
          tcounts1[[i]][[j]] <- Matrix(0, ncol = 1, nrow = length(states[[ti]]))
        } else {
          tcounts1[[i]][[j]] <- Matrix(0, ncol = 1, nrow = 1)
        }
      }
    }
  }
}
for ( i in 1:M) {
  if (any(c(1,26,495,520) == 1)) {
    near.i <- getNbGrids(i)
    ti <- which(lookup == i)
    if(sum(ti) != 0) {
      tcounts2[[i]] <- Matrix(0, ncol = length(near.i) 
                              , nrow = length(states[[ti]]))
    } else {
      tcounts2[[i]] <- Matrix(0, ncol = length(near.i), nrow = 1)
    }
  }
}
print("Init counts...")
########################################################################################
  # Main Loop
  iter <- 0
print("Iter starts...")
  while (iter < numIter) {
    # 1) reset counts
    if (iter > 0) {
      tcounts <- resetCounts1(M, states, lookup)
      tcounts1 <- tcounts[[1]]
      tcounts2 <- tcounts[[2]]
      tcountsg <- tcounts[[3]]
    }
    # Updates counts from paths
    bookL <- lapply(1:M, function(i) lapply(1:length(getNbGrids(i)), function(j) NULL ))
    print("update counts ...")
    for (i in 1:length(paths)) {
      print(i)
      if(length(paths[[i]]$grids) > 5) {
        for (j in 2:length(paths[[i]]$grids)) {
          gn <- paths[[i]]$grids[j]
          ln <- L[[i]][j]
          lc <- L[[i]][j-1]
          gc <- paths[[i]]$grids[j-1]
          gi <- which(lookup == gc)
          gk <- which(getNbGrids(gc) == gn)
          if( any(getNbGrids(gc) == gn) ) {
            tcounts1[[gc]][[gk]][lc, ln] <- tcounts1[[gc]][[gk]][lc, ln] + 1
            if ( class(bookL[[gc]][[gk]]) != "list") {
              if (sum(gi) != 0) {
                if ( (length(states[[gi]]) == 1) && (states[[gi]] == 1) ){ 
                  bookL[[gc]][[gk]] <- lapply(1:length(states[[gi]]), function(i) 1)
                } else {
                  bookL[[gc]][[gk]] <- lapply(1:length(states[[gi]]), function(i) NULL)
                }
              } else {
                bookL[[gc]][[gk]] <- lapply(1:1, function(i) NULL)
              }
            } else {
              bookL[[gc]][[gk]][[lc]] <- c(bookL[[gc]][[gk]][[lc]], ln)
            }
            tcounts2[[gc]][lc, gk] <- tcounts2[[gc]][lc, gk] + 1
          } else {
            tcountsg[gc, gn] <- tcountsg[gc, gn] + 1
          }
        }
      }
    }
    tcounts <- list(tcounts1, tcounts2, tcountsg)
    # Forward Filtering and Backward Sampling
    print("pFFBS starts..")
    sfExport("tcounts", "L","new1000")
    L <- sfLapply(1:N, pFFBSIhmm, paths, tcounts, L, new1000, lookup)
    sfRemove("tcounts", "L","new1000")
    print("pFFBS ends..")
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
    new1000 <- list()
    new1000$CK <- states
    new1000$phi <- phi
    new1000$d <- d
    new1000$beta.k <- bk
    new1000$lpi <- lpi
    new1000$L <- Lhdp
    iter <- iter + 1
  }
  saveRDS(tcounts,"results/tcounts20.rds")
  saveRDS(L,"results/Lf20.rds")
  saveRDS(new1000, "results/new10020.rds")
  sfStop()

