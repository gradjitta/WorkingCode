resultsScriptn <- function(M, sSize,sample, steps, paths, position, startcount, ptrans1, cSize, Lmin, Bmin) {
   ####
   # corridor###
   # For 8
   pat8 <- c( 1, 2,  3, 4,  5)#,10,13,13,131,34,475)
   pos8 <- c(3, 3, 3, 3, 3)#,74, 8,28,  5,12, 27)
   #####
   pos <- pos8
   pat <- pat8
   L0 <- length(pat)
   resl <- matrix(0, L0, sSize)
   resm <- matrix(0, L0, sSize)
   for (j in 1:length(pat)) {
      res  <- lapply(1:sSize, function(i) predHMMdtc(M,sample, steps, pat[j], pos[j], startcount, ptrans1, cSize, Lmin, Bmin))
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
   ret
}
