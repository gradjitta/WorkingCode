getTransg <- function(g, pcounts1, pcounts2, pcountsg, L, K0, lpi, lookup) {
  #getTrans(paths[[i]]$grids, tcounts1, tcounts2, tcountsg,
  #                             L[[i]], sample.hdp$CK, sample.hdp$lpi, lookup)
  set.seed(46)
  N <- length(g)
  MC <- 1e-100
  tempMat3 <- (as.matrix(pcountsg) + 1)
  pcountsg <- t((apply(tempMat3,1,function(x) colMeans(rdirichlet(10,x)))))
  tr.probs <- rep(0,(N-1))
  transMat <- lapply(1:(N-1), function(i) NULL)
  u <- rep(0, N-1)
  for(j in 2:N) {
    print(j)
    time2 <- Sys.time()
    gn <- g[j]
    ln <- L[j]
    lc <- L[j-1]
    gc <- g[j-1]
    if( any(getNbGrids(gc) == gn) ) {
      #t1.pmat <- matrix(0, ncol=t2, nrow=t1)
      #t2.pmat <- matrix(0, ncol=length(getNbGrids(gc)), nrow=t1)
      gk <- which(getNbGrids(gc) == gn)
      K1 <-dim(pcounts1[[gc]][[gk]])[1]
      tempMat1 <- (as.matrix(pcounts1[[gc]][[gk]]) )
      t1.pmat <- t((apply(tempMat1,1,function(x) rdirichlet(1,(x+lpi[[gc]]))+ MC )))
      K2 <-dim(pcounts2[[gc]])[1]
      tempMat2 <- (as.matrix(pcounts2[[gc]]) )
      t2.pmat <- t((apply(tempMat2,1,function(x) rdirichlet(1,(x+lpi[[gc]]))+ MC )))
      #t1.pmat <- pcounts1[[gc]][[gk]]
      #t2.pmat <- pcounts2[[gc]]
      gk <- which(getNbGrids(gc) == gn)
      #tr.probs[j-1] <- t1.pmat[lc, ln] * t2.pmat[lc, gk]
      temp.mat <- matrix(0, ncol=K0, nrow= K0)
      # Slow version
      #for(kk in 1:t2) {
      #  temp.mat[ , kk] <- t1.pmat[, kk] * t2.pmat[, gk]
      #}
      temp.mat <- t1.pmat * t2.pmat[, gk]
      #temp.mat <- t(sweep(t(temp.mat), 2, rowSums(temp.mat), "/"))
      transMat[[j-1]] <- temp.mat / rowSums(temp.mat)
      tr.probs[j-1] <- temp.mat[lc, ln]
      #print(Sys.time() - time2)
    } else {
      tr.probsg <- pcountsg[gc, ]
      sv <- which(lookup == gn)
      if (sum(sv)==0) {
        temp <- rep(1e-10, K0)
      } else {
        temp <- lpi[[sv]]
      }
      tr.probs[j-1] <- (tr.probsg[gn] * temp[ln]) / length(temp)
      # - - - - - - - - - - - - - -
      temp.mat <- matrix(0, ncol=K0, nrow= K0)
      for(kk in 1:K0) {
        temp.mat[ , kk] <- rep((tr.probsg[gn] * temp[kk]) , K0)
      }
      #  - - - - - one - - - - 
      #temp.mat <- t(matrix(1, ncol=t1, nrow= t2) * (temp*tr.probsg[gn]))
      #  - - - - - two - - - - 
      #temp.mat <- matrix(rep((tr.probsg[gn] * temp), rep(t1, t2)), ncol = t2, nrow = t1)
      #temp.mat <- t(sweep(t(temp.mat), 2, rowSums(temp.mat), "/"))
      #transMat[[j-1]] <- temp.mat
      #t(matrix(1, ncol=t1, nrow= t2) * (temp*tr.probsg[gn]))
      # - - - - - - - - - - - 
      transMat[[j-1]] <- temp.mat / rowSums(temp.mat)
    }
    u[j-1] <- runif(1, 0, tr.probs[j-1])
    #print(paste(t1, t2))
    print(Sys.time() - time2)
  }
  ret <- list()
  ret$trans <- tr.probs
  ret$transMat <- transMat
  ret$u <- u
  return(ret)
}