simpleModel <- function(M, sample1, steps, startpos, pathid) {
  library("gtools")
  library("Matrix")
  utcounts <- sample1$counts
  predictDts <- function(steps, jumpcounts, Ntcounts) {
   posi <- pos
   temp <- 0
   gall <- c()
   while(temp < steps) {
     flag <- rbinom(1,size=1,prob = jumpProbs[gt])
     gt <- predNexts(gt, flag, Ntcounts, jumpcounts)
     gall <- c(gall, gt)
     temp <- temp + 1
     #print(exp(tsample))
   }
   #indx <- (pos):length(paths1[[idp]]$grids)
   #tidx <- which(cumsum(exp(sapply(1:length(indx), function(i) paths1[[idp]]$emit[[indx[i]]][1]))) < dt)
   print("pred grid states")
   print(paste(gall))
   print(paste("true states ",steps,"steps ahead"))
   print(paste(paths1[[idp]]$grids[(pos+1):(pos+steps)]))
   ret <- gall[length(gall)]
   ret
 }
 predNexts <- function(gt, flag, Ntcounts, jumpcounts) {
   if (flag == 0) {
     nextgi <- sample(1:length(Ntcounts[[gt]]), size = 1, prob = rdirichlet(1, Ntcounts[[gt]]+1 ))
     gn <- getNbGrids(gt)[nextgi]
   } else {
     gindx <- 1:M
     nextgi <- sample(1:length(jumpcounts[gt, ]), size = 1, prob = rdirichlet(1, jumpcounts[gt, ]+1))
     gn <- gindx[nextgi]
   }
   gn
 }
 Ntcounts <- lapply(1:M, function(i) colSums(utcounts[[2]][[i]]))
 idp <- pathid
 dt <- steps
 pos <- startpos
 Jumps <-  rowSums(utcounts[[3]])
 nonJumps <- sapply(1:M, function(i) sum(colSums(utcounts[[2]][[i]])))
 jumpProbs <- mapply(function(i,j) rbeta(1,1+i,1+j), Jumps, nonJumps)
 gt <- paths1[[idp]]$grids[pos]
 ret <-  predictDts(dt, utcounts[[3]], Ntcounts)
 ret
}
