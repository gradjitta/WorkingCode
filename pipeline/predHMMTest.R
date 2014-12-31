predHMMTest <- function(M, sample, stepsahead, pathid, startpos, startcount, ptrans, cSize, Lmin, Bmin) {
 source("alphaInitgc.R")
 source("getTransgnpTest.R")
 source("simpleModel.R")
 source("getNbGridsF.R")
 source("getNbGrids.R")
 source("getEmitgt.R")
 source("getDist.R")
 source("getBoundingBoxc.R")
 resinit <- sample
 N1 <- ncol(resinit$phi) / 2
 idp <- pathid
 dt <- stepsahead
 pos <- startpos# resinit <- readRDS("../../../results/resultsalpha3g100.rds")
 Jumps <-  rowSums(resinit$counts[[3]])
 nonJumps <- sapply(1:M, function(i) sum(colSums(resinit$counts[[2]][[i]])))
 jumpProbs <- mapply(function(i,j) rbeta(1,1+i,1+j), Jumps, nonJumps)
 alphaTrans <- 10 
 L <- resinit$L
 trans.list <- getTransgnpTest(paths1[[idp]]$grids, resinit$probs[[1]], resinit$probs[[2]], 
                           resinit$probs[[3]], L[[idp]], resinit$K-1, resinit$lw, lookup1, 
                           resinit$K, alphaTrans, cSize)
 trans.path <- trans.list$transMat
 emit.path <- getEmitgt(resinit$K, paths1[[idp]]$emit, paths1[[idp]]$grids, resinit$phi, lookup1, N1)
 alphas <- alphaInitgc(paths1[[idp]], L[[idp]], lookup1, resinit$phi, resinit$K)
 alpha0 <- startcount / sum(startcount)
 alphatemp <- alpha0[paths1[[idp]]$grids[1]] * resinit$lw[paths1[[idp]]$grids[1], ]
 alpha0 <- alphatemp / sum(alphatemp)
 getAlphas <- function(trans.path, emit.path, grids, u, T0) {
     for (t in 2:T0) {
       temp.trans <- trans.path[[t-1]]
       #sum.trans <- colSums(sweep((u[t-1] < temp.trans),MARGIN=1,alphas[[t-1]],`*`))
       sum.trans <- colSums(temp.trans * alphas[[t-1]])
       temp.a <- (emit.path[[t]] * sum.trans)
       alphas[[t]] <- temp.a / sum(temp.a)
     }
     alphas
 }
 alphas[[1]] <- alpha0
 alphas <- getAlphas(trans.path, emit.path, paths1[[idp]]$grids, u, length(paths1[[idp]]$grids))

 
# This is the main prediction function, that generates grid blocks and clusters used given previous.
predNextii<- function(ptrans,phi, alphas, gt,lt, lookup1,K, lw,t,dt, flag, cSize, Lmin, Bmin) {
   gn <- getNbGridsF(gt, cSize, Lmin, Bmin)
   Ln <- K
   if (flag == 0) {    ## No Jumps
     retval <- list()
     # sample next grid block
     sampg <- ptrans[[2]][[gt]][lt, ]
     sampg <- sampg / sum(sampg)
     gi <- sample(1:length(sampg),1,prob = sampg)
     # sample next local cluster used
     sampl <- colSums((ptrans[[1]][[gt]][[gi]] * (ptrans[[2]][[gt]][, gi] * alphas) ))
     alpha0 <- sampl / sum(sampl)
     gl <- sample(1:length(alpha0),1,prob = alpha0)
     retval$g <- gn[gi]
     retval$l <- gl
     retval$alphas <- alpha0
   } else {    ## Jumps
     retval <- list()
     lnext <- sample(1:length(lw[gt, ]), 1, prob = lw[gt, ])
     jprobs <- ptrans[[3]][gt, ] * lw[gt, lnext]
     gnext <- sample(1:length(jprobs), 1,  prob = jprobs)
     retval$g <- gnext
     retval$l <- lnext
     retval$alphas <- ptrans[[3]][gt, gnext] * lw[gt, ] / sum(ptrans[[3]][gt, gnext] * lw[gt, ])
  }
   return (retval)
}


 predictDt <- function(dt, ptrans, alphas, gt,lt, lookup1,K, lw,t,phi, cSize, Lmin, Bmin) {
 # This function calls predNextii() 
   temp <- 0
   gall <- c()
   lall <- c()
   while(temp < dt) {
     flag <- rbinom(1,size=1,prob = jumpProbs[gt])
     retvals <- predNextii(ptrans,phi, alphas, gt,lt, lookup1,K, lw,t,dt,flag, cSize, Lmin, Bmin)
     gt <- retvals$g
     gall <- c(gall, gt)
     lt <- retvals$l
     lall <- c(lall, lt)
     alphas <- retvals$alphas
     temp <- temp + 1
   }
   gtr <- paths1[[idp]]$grids[pos + dt]
   ret <- c(gall[length(gall)], gtr, length(gall))
   ret
 }
 Ltemp <- sample(1:length(alphas[[pos]]),1,prob = alphas[[pos]])
 gret1 <- predictDt(dt, ptrans, alphas[[pos]], paths1[[idp]]$grids[pos], Ltemp, 
                    lookup1, resinit$K, resinit$lw, pos, resinit$phi, cSize, Lmin, Bmin)
 # simpleModel:: First order Markov Model.
 gret2 <- simpleModel(M, resinit, dt, pos, idp, dt, cSize, Lmin, Bmin)
 dists <-  getDist(gret1[1], gret2, gret1[2], cSize)
 print(dists)
}
