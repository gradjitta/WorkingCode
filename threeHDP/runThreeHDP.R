library(Matrix)
library(parallel)
library("gtools")
source("resetCountsK.R")
source("initStatesSimp.R")
source("dataForfhmm.R")
source("sampleMeanjnm.R")
source("sampleVariancejn.R")
source("updateCountsS.R")
source("getTransProbsfc.R")
source("pFFBSIhmmfnp.R")
source("getTransgnp.R")
source("getEmitg.R")
source("alphaInitg.R")
source("getNbGrids.R")
source("getMCounts.R")
source("getUsedClusters.R")
source("sampleTableCount.R")
source("stirlingNumbers.R")
source("updateGlobalWeights.R")
source("updateCellWeights.R")
lookup1 <- readRDS("newData/lookup1.rds")
paths1 <- readRDS("newData/paths.rds")
new1000 <- readRDS("newData/new1000c.rds")
gridlist <- readRDS("newData/gridlist.rds")
print("create stirling numbers")
sN <- stirlingNumbers(1000)
# parameters/ hyperparameters
print("init parameters")
mean.init <- 1.5
s <- 1
gamma <- 5
alpha <- 10
alphaTrans <- 10
K <- 40
nCell <- 520
N0 <- 5
print("init L, counts, nK and m")
globalWeights <- rdirichlet(1,rep(gamma,K))
cellWeights <- matrix(0,nCell,K)
L <- initStatesSimp(lookup1, paths1, K)
tcountsi <- resetCountsK(nCell, K, lookup1)
utcounts <- updateCountsS(nCell, gridlist, paths1, lookup1, K, L, tcountsi)
nK <- getMCounts(nCell, utcounts[[1]], K)
m <- matrix(0,nCell,K)
for(gt in 1:nCell) {
  nTemp <- nK[[gt]]
  cellWeights[gt,] <- rdirichlet(1,alpha*globalWeights)
  out <- updateCellWeights(nTemp,cellWeights[gt,],globalWeights,alphaTrans,alpha,sN)
  cellWeights[gt,] <- out$cellWeights
  m[gt,] <- out$m
}
print("##################")
nIter <- 1
for(iter in 1:nIter) {
  print(paste(" iteration",iter))
  #
  # Update the global weights
  #
  print("update global weights..")
  globalWeights <- updateGlobalWeights(m,globalWeights,alpha,gamma,sN,maxNew=2)
  Kold <- K
  K <- length(globalWeights)
  print(globalWeights)
  # Sample parameters for the newly created sticks
  # Update Phi
  print("data for phi")
  data.hmm <- dataForfhmm(paths1, L, lookup1)
  d <- data.hmm$Lhdp
  data <- data.hmm$data
  J <- dim(summary(data))[1]
  print("update phi")
  phi <- matrix(0, K, 2*N0)
  for (j in 1:K) {
    ###################
    if(j<=Kold) {
      # sample mean with prior Normal(0, 1/s) on mean
      phi[j, 1:N0] <- sampleMeanjnm(data, d, j,
                                  phi[j,(N0+1):(2*N0),drop = T], s ,mean.init, J)
      # sample 1/variance with prior Gamma(epi, epi) on 1/variance
      phi[j, (N0+1):(2*N0)] <- sampleVariancejn(data,
                                              j, 1, d, phi[j, 1:N0, drop = T], J)
    } else {
      m0 <- rnorm(N0, mean.init, sqrt(1/s))
      v0 <- 1/rgamma(N0, 1, 1)
      phi[j, ] <- c(m0, v0)
    }
  }
  #
  # Update the other sets of weights
  #
  print("update cell weights..")
  m <- matrix(0,nCell,K)
  cellWeightsOld <- cellWeights
  cellWeights <- matrix(0,nCell,K)
  for(gt in 1:nCell) {
    # Update the target cell weights and the numbers of tables
    # at the lower level
    
    # Fancy way of creating a matrix that contains all transitions
    # into this cell; it does not matter in which order they are
    nTemp <- nK[[gt]]
    out <- updateCellWeights(nTemp,cellWeightsOld[gt,],globalWeights,alphaTrans,alpha,sN)
    cellWeights[gt, ] <- out$cellWeights
    m[gt, ] <- out$m
  }
  
  # Update assignments using FB and then re-compute the counts
  # based on the likelihood and the transition matrix
  # pFFBSIhmmfnp(100, paths1, trprobs, L, new1000, lookup1, lpi42, 40, phi, 42, 10)
  print("get probs..")
  pcounts <- getTransProbsfc(utcounts, cellWeights, alphaTrans, K)
  print("FFBS..")
  L <- mclapply(1:length(paths1), pFFBSIhmmfnp, paths1, pcounts, L, lookup1, cellWeights, Kold, phi, K, alphaTrans, mc.cores = 8)
  print("update Counts.. ")
  tcountsi <- resetCountsK(nCell, K, lookup1)
  utcounts <- updateCountsS(nCell, gridlist, paths1, lookup1, K, L, tcountsi)
  nK <- getMCounts(nCell, utcounts[[1]], K)
  # Kill unnecessary components
  print("throw unnecessary components..")
  if(iter<nIter) {
    keep <- which(getUsedClusters(nCell, ptcounts[[1]], K) > 0)
    cellWeights <- cellWeights[,keep]
    globalWeights <- globalWeights[keep]
    phi <- phi[keep, ]
    m <- m[keep, ]
    K <- length(keep)
  }
}
