trans.list <- getTransp(paths[[i]]$grids, tcounts1, tcounts2, tcountsg, tempL, sample.hdp$CK, sample.hdp$lpi, lookup)
trans.path <- trans.list$transMat
emit.path <- getEmitf(paths[[i]]$emit, paths[[i]]$grids, tempL,sample.hdp$CK, sample.hdp$phi, lookup)
alphas <- alphaInit(paths[[i]], L20, new1000c$CK, lookup1, new1000c$phi)

getAlphas <- function(trans.path, emit.path, grids, u) {
    for (t in 2:50) {
      temp.trans <- trans.path[[t-1]]
      #sum.trans <- colSums(sweep((u[t-1] < temp.trans),MARGIN=1,alphas[[t-1]],`*`))
      sum.trans <- colSums( ((u[t-1] < temp.trans) * alphas[[t-1]]) )
      if(sum(emit.path[[t]]) == 0) {
        temp.a <- (sum.trans)
      } else {
        temp.a <- (emit.path[[t]] * sum.trans)
      }
      if (sum(temp.a) == 0) {
        temp.a <- rep(1, length(sum.trans))
      }
      alphas[[t]] <- temp.a / sum(temp.a)
    }
    if(sum(emit.path[[1]]) != 0) {
      alphas[[1]] <- emit.path[[1]]
    }
    alphas
}



predNext <- function(ptrans, alphas, gt,lt, new1000c, lookup, N) {
  gn <- getNbGrids(gt)
  Ln <- length(new1000c$CK[[which(gt == lookup)]])
  nextVals <- lapply(1:length(gn), function(i) NULL)
  for (i in 1:length(gn)) {
   nextVals[[i]] <- ptrans[[1]][[gt]][[i]][lt, ] * ptrans[[2]][[gt]][lt, i] * alphas[[N]][lt]
  }
  maxvalind <- sapply(1:length(nextVals), function(j) which(nextVals[[j]] == max(nextVals[[j]])) )
  maxvals <- sapply(1:length(nextVals), function(j) max(nextVals[[j]]) )
  lnext <- maxvalind[which(maxvals == max(maxvals))]
  gnext <- gn[which(maxvals == max(maxvals))]
  print(paste(gnext, lnext))
}


