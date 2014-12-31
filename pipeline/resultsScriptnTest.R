resultsScriptnTest <- function(M, sSize,sample, steps, paths, position, startcount, ptrans1,ptrans2, cSize, Lmin, Bmin, usePaths) {
   ####
   # corridor###
   # For 8
   pat <- 201:300
   pos <- rep(5,100)
   #####
   #pos <- 32
   #pat <- 1
   # For 30
   pos <- rep(3,100)
   pat <- usePaths[1:100]
   reslist <- list()
   L0 <- length(pat)
   resl <- matrix(0, L0, sSize)
   resm <- matrix(0, L0, sSize)
   posSet <- c(3)
   #posSet <- c(2,3,4,6)
   for (ii in 1:length(posSet)) {
   nstart <- posSet[ii]
   pos <- rep(nstart, 100)
   for (j in 1:length(pat)) {
      res  <- lapply(1:sSize, function(i) predHMMTest(M,sample, steps, pat[j], pos[j], startcount, ptrans1,ptrans2, cSize, Lmin, Bmin))
      localmodel <- sapply(1:length(res), function(i) res[[i]][1])
      x1 <- sapply(1:sSize,function(i) sum(localmodel <= (i-1) ))
      markovmodel <- sapply(1:length(res), function(i) res[[i]][2])
      x2 <- sapply(1:sSize,function(i) sum(markovmodel <= (i-1) ))
      resl[j, ] <- x1
      resm[j, ] <- x2
   }
   ret <- list()
   ret$th <- 1:sSize
   ret$lm <- resl
   ret$mm <- resm
   reslist[[ii]]<-ret
   }
   reslist
}
