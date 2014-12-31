simpleModel <- function(M, sample1, timeahead, startpos, pathid, steps, cSize, Lmin, Bmin) {
  library("gtools")
  library("Matrix")
  #paths1 <- readRDS("newData/lpathsb.rds")
  #baskets <- readRDS("newData/basketsn.rds")
  #lookup1 <- baskets$lookup1
  #gridlist <- baskets$gridlist
  utcounts <- sample1$counts
  predictDts <- function(dt, jumpcounts, Ntcounts) {
   posi <- pos
   temp <- 0
   gall <- c()
   while(temp < steps) {
#     print(gt)
     flag <- rbinom(1,size=1,prob = jumpProbs[gt])
     gt <- predNexts(gt, flag, Ntcounts, jumpcounts)
     gall <- c(gall, gt)
     temp <- temp + 1
     #print(exp(tsample))
   }
   indx <- (pos):length(paths1[[idp]]$grids)
  # tidx <- which(cumsum(exp(sapply(1:length(indx), function(i) paths1[[idp]]$emit[[indx[i]]][1]))) < dt)
   print("pred grid states")
   print(paste(gall))
   print(paste("true states ",dt,"seconds ahead"))
#   print(paste(paths1[[idp]]$grids[pos + steps -1]))
   print(paste(paths1[[idp]]$grids[pos + steps]))
   ret <- gall[length(gall)]
   ret
 }
 predNexts <- function(gt, flag, Ntcounts, jumpcounts) {
   if (flag == 0) {
     
     nextgi <- sample(1:length(Ntcounts[[gt]]), size = 1, prob = rdirichlet(1, Ntcounts[[gt]]+1 ))
     gn <- getNbGridsF(gt, cSize, Lmin, Bmin)[nextgi]
   } else {
     gindx <- 1:M
     nextgi <- sample(1:length(jumpcounts[gt, ]), size = 1, prob = rdirichlet(1, jumpcounts[gt, ]+1))
     gn <- gindx[nextgi]
   }
   gn
 }
 Ntcounts <- lapply(1:M, function(i) colSums(utcounts[[2]][[i]]))
 idp <- pathid
 dt <- timeahead
 pos <- startpos
 Jumps <-  rowSums(utcounts[[3]])
 nonJumps <- sapply(1:M, function(i) sum(colSums(utcounts[[2]][[i]])))
 jumpProbs <- mapply(function(i,j) rbeta(1,1+i,1+j), Jumps, nonJumps)
 gt <- paths1[[idp]]$grids[pos]
 ret <-  predictDts(dt, utcounts[[3]], Ntcounts)
 ret
}
