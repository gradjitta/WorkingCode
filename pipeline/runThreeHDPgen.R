set.seed(42)
library("Matrix")
library("parallel")
library("gtools")
source("resetCountsKc.R")
source("initStatesSimp.R")
source("dataForfhmmt.R")
source("sampleMeanjnm.R")
source("sampleVariancejn.R")
source("updateCountsS.R")
source("getTransProbsfcp.R")
source("pFFBSIhmmfnpt.R")
source("getTransgnp.R")
source("getEmitgt.R")
source("alphaInitgc.R")
source("getNbGridsF.R")
source("getMCounts.R")
source("getUsedClusters.R")
source("sampleTableCount.R")
source("stirlingNumbers.R")
source("updateGlobalWeights.R")
source("updateCellWeights.R")
source("getMCountsc.R")
source("updateCountsS.R")
source("updateL.R")
source("getCorners.R")
##################
baskets <- readRDS("data/paths30pmcData.rds")
lookup1 <- baskets[[5]]
paths1 <- baskets[[2]]
gridlist <- baskets[[4]]
################
#lookup1 <- readRDS("newData/lookup1.rds")
#paths1 <- readRDS("newData/paths.rds")
#new1000 <- readRDS("newData/new1000c.rds")
#gridlist <- readRDS("newData/gridlist.rds")
print("create stirling numbers")
sN <- stirlingNumbers(9000)
# parameters/ hyperparameters
print("init parameters")
mean.init <- 1.5
s <- 1
#gamma <- 5
#alpha <- 10
#alphaTrans <- 10
Lmin <- 100
Bmin <- 80
cSize <- 30
nCell <- (floor(Bmin / cSize) +1) * (floor(Lmin / cSize) + 1)
gamma <- 10
alpha <- 3
alphaTrans <- 10
K <- 1
keep <- 1:K
N0 <- 12
print("init L, counts, nK and m")
globalWeights <- rdirichlet(1,rep(gamma,K))
cellWeights <- matrix(0,nCell,K)
# L <- initStatesSimp(lookup1, paths1, K)
phi <- matrix(0, K, 2*N0)
for (i in 1:K) {
  m0 <- rnorm(N0, mean.init, sqrt(1/s))
  v0 <- 1/rgamma(N0, 1, 1)
  phi[i, ] <- c(m0, v0)
}
# Init with kmeans
#L <- readRDS("newData/Lkmeans.rds")
#L <- readRDS("newData/L5.rds")
#L <- readRDS("newData/L1.rds")
L <- lapply(1:length(paths1), function(i) rep(1,length(paths1[[i]]$grids)))
#L <- readRDS("newData/L10.rds")
#tcountsi <- resetCountsKc(nCell, K, lookup1)
tcountsi <- resetCountsKc(nCell, K, lookup1, cSize, Lmin, Bmin)
utcounts <- updateCountsS(nCell, gridlist, paths1, lookup1, K, L, tcountsi,cSize)
# NEED TO CHANGE getMCounts to use jumpCounts
nK <- getMCounts(nCell, utcounts[[1]], keep,cSize)
m <- matrix(0,nCell,K)
for(gt in 1:nCell) {
  nTemp <- nK[[gt]]
  cellWeights[gt,] <- rdirichlet(1,alpha*globalWeights)
  out <- updateCellWeights(nTemp,cellWeights[gt,],globalWeights,alphaTrans,alpha,sN)
  cellWeights[gt,] <- out$cellWeights
  m[gt,] <- out$m
}
print("##################")
nIter <- 1200
for(iter in 1:nIter) {
  print(paste(" iteration",iter))
  #
  # Update the global weights
  #
  if(iter == 1){
    print(getUsedClusters(nCell, utcounts[[1]], K) + utcounts[[4]])
  }
  print("update global weights..")
  globalWeights <- updateGlobalWeights(m,globalWeights,alpha,gamma,sN,maxNew=1)
  Kold <- K
  K <- length(globalWeights)
  print(globalWeights)
  # Sample parameters for the newly created sticks
  # Update Phi
  print("data for phi")
  data.hmm <- dataForfhmmt(paths1, L, lookup1, N0, nCell)
  d <- data.hmm$Lhdp
  data <- data.hmm$data
  J <- dim(summary(data))[1]
  print("update phi")
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
      phi <- rbind(phi,  c(m0, v0))
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
  pcounts <- getTransProbsfcp(nCell, utcounts, cellWeights, alphaTrans, K, cSize)
  print("FFBS..")
  alphacounts <- sapply(1:K, function(j) sum( sapply(1:length(L), function(i) L[[i]][1]) == j))
  alpha0 <- rdirichlet(1, alphacounts+1)
  alpha0 <- alpha0[,,drop = T]
  startcount <- sapply(1:nCell, function(i) sum(sapply(1:length(paths1), function(j) paths1[[j]]$grids[1] == i))) 
  alpha0 <- startcount / sum(startcount)
  L <- mclapply(1:length(paths1), pFFBSIhmmfnpt, paths1, pcounts, L, lookup1, cellWeights, Kold, phi, K, alphaTrans, alpha0, cSize,mc.cores = 8)
  print("update Counts.. ")
  tcountsi <- resetCountsKc(nCell, K, lookup1, cSize, Lmin, Bmin)
  utcounts <- updateCountsS(nCell, gridlist, paths1, lookup1, K, L, tcountsi, cSize)
  alphacounts <- sapply(1:K, function(j) sum( sapply(1:length(L), function(i) L[[i]][1]) == j))
  # Kill unnecessary components
  print("save every 10th")
  if (iter %% 1 == 0) {
    results <- list()
    results$iter <- iter
    results$K <- K
    results$L <- L
    results$m <- m
    results$n <- nK
    results$phi <- phi
    results$gw <- globalWeights
    results$lw <- cellWeights
    results$probs <- pcounts
    results$counts <- utcounts
    saveRDS(results,"data/results/resultsalpha3g10C30Trainpmcl.rds")
  }
  print("throw unnecessary components..")
  if(iter<nIter) {
    ckeep <- getUsedClusters(nCell, utcounts[[1]], K) + utcounts[[4]] + alphacounts
    print(ckeep)
    keep <- which(ckeep > 0)
    cellWeights <- cellWeights[,keep, drop = F]
    globalWeights <- globalWeights[,keep, drop = F]
    phi <- phi[keep, ,drop = F]
    m <- m[, keep, drop = F]
    L <- updateL(L, keep)
    if (K != length(keep)) {
      K <- length(keep)
      tcountsi <- resetCountsKc(nCell, K, lookup1,cSize, Lmin, Bmin)
      utcounts <- updateCountsS(nCell, gridlist, paths1, lookup1, K, L, tcountsi, cSize)
      nK <- getMCountsc(nCell, utcounts[[1]], keep, cSize)
    } else {
      nK <- getMCountsc(nCell, utcounts[[1]], keep, cSize)
    }
    K <- length(keep)
  }
}
results <- list()
results$K <- K
results$L <- L
results$m <- m
results$n <- nK
results$phi <- phi
results$gw <- globalWeights
results$lw <- cellWeights
results$probs <- pcounts
results$counts <- utcounts
saveRDS(results,"data/results/resultsThree800alpha3g100C30Trainpmcl.rds")
