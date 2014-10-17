source("alphaInitgc.R")
source("getTransgnp.R")
library("gtools")
resinit <- sample
idp <- pathid
timeahead <- dt
startpos <- pos
paths1 <- readRDS("newData/lpathsb.rds")
baskets <- readRDS("newData/basketsn.rds")
lookup1 <- baskets$lookup1
gridlist <- baskets$gridlist
Jumps <-  rowSums(resinit$counts[[3]])
nonJumps <- sapply(1:520, function(i) sum(colSums(resinit$counts[[2]][[i]])))
jumpProbs <- mapply(function(i,j) rbeta(1,1+i,1+j), Jumps, nonJumps)
endPos <- sapply((1+26):(26+26), function(i) sum(sapply(1:length(resinit$L),function(j) paths1[[j]]$grids[length(resinit$L[[j]])] == i)) )
endindex <- c(27:52)
endprobs <- rep(0,520)
endprobs[endindex] <- endPos
endprobs <- endprobs / sum(endprobs)
alphaTrans <-10
L <- resinit$L
trans.list <- getTransgnp(paths1[[idp]]$grids, resinit$probs[[1]], resinit$probs[[2]], resinit$probs[[3]], L[[1]], resinit$K-1, resinit$lw, lookup1, resinit$K, alphaTrans)
trans.path <- trans.list$transMat
emit.path <- getEmitg(resinit$K, paths1[[idp]]$emit, paths1[[idp]]$grids, resinit$phi, lookup1)
L <- resinit$L
u <- trans.list$u
alphacounts <- sapply(1:resinit$K, function(j) sum( sapply(1:length(L), function(i) L[[i]][1]) == j))
alpha0 <- rdirichlet(1, alphacounts+1)
alpha0 <- alpha0[,,drop = T]
alphas <- alphaInitgc(paths1[[idp]], L[[idp]], lookup1, resinit$phi, resinit$K)
alphas[[1]] <- alpha0
getAlphas <- function(trans.path, emit.path, grids, u, T0) {
    for (t in 2:T0) {
      temp.trans <- trans.path[[t-1]]
      #sum.trans <- colSums(sweep((u[t-1] < temp.trans),MARGIN=1,alphas[[t-1]],`*`))
      sum.trans <- colSums( ((u[t-1] < temp.trans) * alphas[[t-1]]) )
      temp.a <- (emit.path[[t]] * sum.trans)
      alphas[[t]] <- temp.a / sum(temp.a)
    }
    alphas
}
alphas <- getAlphas(trans.path, emit.path, paths1[[1]]$grids, u, length(paths1[[1]]$grids))
getIndMat <- function(tind, tlen, nr, nc) {
   temp <- rep(0, tlen)
   temp[tind] <- 1
   tempmat <- matrix(temp, nrow = nr, ncol = nc)
   ret <- which(tempmat == 1, arr.ind = TRUE)
   ret
}

getJumpMat <- function(gt, probs, lw, phi, alphas) {
   K <- dim(lw)[2]
   M <- dim(lw)[1]
   jumpMat <- matrix(0, K, M)
   for (i in 1:K) {
    tempemit <- rnorm(5,phi[i,1:5], sqrt(phi[i,6:10]))
    emits <- tempEmit(K, tempemit, phi)
    jumpMat[i, ] <- probs[gt, ] * lw[, i]
   }
   jumpMat <- jumpMat * alphas
}
tempEmit <- function(K, emit, phi) {
  N <- dim(phi)[2] / 2
  m <- 1:K
  temp.emit <- rep(0, length(m))
  for(j in 1:length(m)) {
    mj <- m[j]
    temp.emit[j] <- sum(dnorm(emit, phi[mj, 1:5], sqrt(phi[mj, 6:10]), log = TRUE))
  }
  num.temp <- exp(temp.emit - max(temp.emit))
  emit.probs <- num.temp / sum(num.temp)
  return(emit.probs)
}
predNext <- function(ptrans,phi, alphas, gt,lt, lookup1,K, lw,t,dt, flag) {
  gn <- getNbGrids(gt)
  Ln <- K
  nextVals <- matrix(0, nrow = length(gn), ncol = Ln )
  for (i in 1:length(gn)) {
   nextVals[i, ] <- colSums(((ptrans[[1]][[gt]][[i]] * ptrans[[2]][[gt]][, i])) * alphas)
  }
  for (l in 1:Ln) {
    tempemit <- rnorm(5,phi[l,1:5], sqrt(phi[l,6:10]))
    emits <- tempEmit(K, tempemit, phi)
    nextVals[, l] <- nextVals[, l]*emits[l]
  }
  localvec <- as.vector(nextVals)
  ##########
  #jumpmat <- getJumpMat(gt, ptrans[[3]], lw, phi)
  #jumpvec <- as.vector(jumpmat)
  ##########
  tempvec <- c(localvec)#,jumpvec)
  tempprobs <- tempvec / sum(tempvec)
  K0 <- length(localvec)
  sampleind <-  1 + sum(runif(length(tempprobs), 0, 1) > cumsum(tempprobs))
  if (flag == 0) {
    nextVs <- getIndMat(sampleind, length(tempvec), length(gn), Ln) 
    gnext <- gn[nextVs[1]]
    lnext <- nextVs[2]
    retval <- list()
    k <- which(gn == gnext)
    retval$g <- gnext
    retval$l <- lnext
    retval$alphas <- nextVals[k, ] / sum(nextVals[k, ])
    print("nj")
  } else {
    print("j")
    retval <- list()
    jumpmat <- getJumpMat(gt, ptrans[[3]], lw, phi, alphas)
    jumpvec <- as.vector(jumpmat)
    jprobs <- jumpvec / sum(jumpvec)
    sampleind <-  1 + sum(runif(length(jprobs), 0, 1) > cumsum(jprobs))
    nextVs <- getIndMat(sampleind, length(jumpvec), Ln, 520)
    gnext <- nextVs[2]
    lnext <- nextVs[1]
    retval$g <- gnext
    retval$l <- lnext
    tempemit <- rnorm(5,phi[lnext,1:5], phi[lnext,6:10])
    emits <- tempEmit(K, tempemit, phi)
    retval$alphas <- jumpmat[, gnext]*emits / sum(jumpmat[, gnext]*emits)
 }
  return (retval)
}
predictDt <- function(dt, ptrans, alphas, gt,lt, lookup1,K, lw,t,phi) {
  temp <- 0
  gall <- c(gt)
  lall <- c(lt)
  while(temp < dt) {
    flag <- rbinom(1,size=1,prob = jumpProbs[gt])
    retvals <- predNext(ptrans,phi, alphas, gt,lt, lookup1,K, lw,t,dt,flag)
    gt <- retvals$g
    gall <- c(gall, gt)
    lt <- retvals$l
    lall <- c(lall, lt)
    alphas <- retvals$alphas
    #temp <- temp + resinit$phi[lt, ][1]
    tsample <- rnorm(1, resinit$phi[lt, ][1],resinit$phi[lt, ][6])
    temp <- temp + exp(tsample)
    #print(exp(tsample))
  }
  tidx <- which(cumsum(exp(sapply(1:length(paths1[[idp]]$grids), function(i) paths1[[idp]]$emit[[i]][1]))) < dt)
  print("pred grid states")
  print(paste(gall))
  print("pred local states")
  print(paste( lall))
  print("true 10 states ahead")
  print(paths1[[idp]]$grids[tidx])
  print(resinit$L[[idp]][tidx])
}
predictDt(dt, resinit$probs, alphas[[pos]], paths1[[idp]]$grids[pos],resinit$L[[idp]][pos] , lookup1,resinit$K, resinit$lw, pos, resinit$phi)
