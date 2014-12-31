library(MCMCpack)
library(Matrix)
library(parallel)
source("resetCountsK.R")
source("initStatesSimp.R")
source("dataForfhmm.R")
source("sampleMeanjnm.R")
source("sampleVariancejn.R")
source("updateCountsS.R")
source("getTransProbsf.R")
source("pFFBSIhmmfin.R")
source("getTransgc.R")
source("getEmitg.R")
source("alphaInitg.R")
source("getNbGrids.R")
lookup1 <- readRDS("newData/lookup1.rds")
paths1 <- readRDS("newData/paths.rds")
new1000 <- readRDS("newData/new1000c.rds")
gridlist <- readRDS("newData/gridlist.rds")
# Init PHI
N0 <- 5
K <- 40
M <- 520
phi <- matrix(0, K, 2*N0)
mean.init <- 1.5
s <- 1
for (i in 1:K) {
  m0 <- rnorm(N0, mean.init, sqrt(1/s))
  v0 <- 1/rgamma(N0, 1, 1)
  phi[i, ] <- c(m0, v0)
}
# init cluster assignments
tcountsi <- resetCountsK(M, K, lookup1)
L <- initStatesSimp(lookup1, paths1, K)
nk <- sapply(1:K, function(j) sum(sapply(1:length(L), function(i) sum(L[[i]] == j) )))
print("init starts")
#nk <- sapply(1:K, function(j) sum(sapply(1:J, function(i) sum(d[[i]] == j))) )
beta.k <- rdirichlet(1, (nk+1))
lpi <- lapply(1:M,function(i) beta.k)
print("data for fhmm")
data.hmm <- dataForfhmm(paths1, L, lookup1)
d <- data.hmm$Lhdp
data <- data.hmm$data
J <- dim(summary(data))[1]
for (j in 1:K) {
    ###################
    # sample mean with prior Normal(0, 1/s) on mean
    phi[j, 1:N0] <- sampleMeanjnm(data, d, j,
                                  phi[j,(N0+1):(2*N0),drop = T], s ,mean.init, J)
    # sample 1/variance with prior Gamma(epi, epi) on 1/variance
    phi[j, (N0+1):(2*N0)] <- sampleVariancejn(data,
                                              j, 1, d, phi[j, 1:N0, drop = T], J)
}
# Main Loop
print("Loop starts")
iter <- 1
totIter <- 5
L <- initStatesSimp(lookup1, paths1, K)
while (iter < totIter) {
    if (iter > 1) {
       print("reset counts..")
       tcountsi <- resetCountsK(M, K, lookup1)
       nk <- sapply(1:K, function(j) sum(sapply(1:length(L), function(i) sum(L[[i]] == j) )))
       #nk <- sapply(1:K, function(j) sum(sapply(1:J, function(i) sum(d[[i]] == j))) )
       beta.k <- rdirichlet(1, (nk+1))
       lpi <- lapply(1:M,function(i) beta.k)
    }
    # update counts
    print("update counts..")
    utcounts <- updateCountsS(M, gridlist, paths1, lookup1, K, L, tcountsi)
    print("get probs..")
    pcounts <- getTransProbsf(utcounts)
    # FFBS
    print("FFBS..")
    time5 <- Sys.time()
    L <- mclapply(1:length(paths1), pFFBSIhmmfin, paths1, pcounts, L, new1000, lookup1, lpi, K, phi, mc.cores = 8)
    print(Sys.time() - time5 )
    # Update d from FFBS
    print("update data..")
    data.hmm <- dataForfhmm(paths1, L, lookup1)
    d <- data.hmm$Lhdp
    data <- data.hmm$data
    J <- dim(summary(data))[1]
    # Update Phi
    for (j in 1:K) {
    ###################
    # sample mean with prior Normal(0, 1/s) on mean
      phi[j, 1:N0] <- sampleMeanjnm(data, d, j,
                                  phi[j,(N0+1):(2*N0),drop = T], s ,mean.init, J)
    # sample 1/variance with prior Gamma(epi, epi) on 1/variance
      phi[j, (N0+1):(2*N0)] <- sampleVariancejn(data,
                                              j, 1, d, phi[j, 1:N0, drop = T], J)
    }
    iter <- iter + 1
}
results <- list()
results$phi <- phi
results$L <- L
results$betak <- beta.k
saveRDS(results, "results/resultsTest.rds")
