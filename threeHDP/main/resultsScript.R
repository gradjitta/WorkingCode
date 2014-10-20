resultsScript <- function(sSize,sample, steps, paths, position, startcount, ptrans1) {
   ####
   # corridor###
   pat <- c( 1, 2, 3, 4 ,5,10,13,13,131,34,475)
   pos <- c(13,13,43,13,24,74, 8,28,  5,12, 27)
   #####
   L0 <- length(pat)
   resl <- matrix(0, L0, sSize)
   resm <- matrix(0, L0, sSize)
   for (j in 1:length(pat)) {
      res  <- lapply(1:sSize, function(i) predHMMdt(sample, steps, pat[j], pos[j], startcount, ptrans1))
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
