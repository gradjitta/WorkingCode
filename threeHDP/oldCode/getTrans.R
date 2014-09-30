getTrans <- function(g, tcounts1, tcounts2, tcountsg, L, Ck, lpi, lookup) {
  #getTrans(paths[[i]]$grids, tcounts1, tcounts2, tcountsg,
  #                             L[[i]], sample.hdp$CK, sample.hdp$lpi, lookup)
  N <- length(g)
  tr.probs <- rep(0,(N-1))
  transMat <- lapply(1:(N-1), function(i) NULL)
  u <- rep(0, N-1)
  for(j in 2:N) {
    gn <- g[j]
    ln <- L[j]
    lc <- L[j-1]
    gc <- g[j-1]
    if (sum(which(lookup == gn)) == 0) {
      tempn2 <- 1
    } else {
      tempn2 <- Ck[[which(lookup == gn)]]
    }
    if (sum(which(lookup == gc)) == 0) {
      tempnc <- 1
    } else {
      tempnc <- Ck[[which(lookup == gc)]]
    }
    t2 <- length(tempn2)
    t1 <- length(tempnc)
    if( any(getNbGrids(gc) == gn) ) {
      t1.pmat <- matrix(0, ncol=t2, nrow=t1)
      t2.pmat <- matrix(0, ncol=length(getNbGrids(gc)), nrow=t1)
      gk <- which(getNbGrids(gc) == gn)
      for(k in 1:t1) {
        t1.pmat[k, ] <- SampleDirichlet(tcounts1[[gc]][[gk]][k, ] + 1)
        t1.pmat[k, ] <- t1.pmat[k, ] / sum(t1.pmat[k, ])
        t2.pmat[k, ] <- SampleDirichlet(tcounts2[[gc]][k, ] + 1)
        t2.pmat[k, ] <- t2.pmat[k, ] / sum(t2.pmat[k, ])
      }
      gk <- which(getNbGrids(gc) == gn)
      #tr.probs[j-1] <- t1.pmat[lc, ln] * t2.pmat[lc, gk]
      temp.mat <- matrix(0, ncol=t2, nrow= t1)
      for(kk in 1:t2) {
        temprow <- t1.pmat[, kk] * t2.pmat[, gk]
        temp.mat[ , kk] <- t1.pmat[, kk] * t2.pmat[, gk]
      }
      temp.mat <- t(sweep(t(temp.mat), 2, rowSums(temp.mat), "/"))
      transMat[[j-1]] <- temp.mat
      tr.probs[j-1] <- temp.mat[lc, ln]
    } else {
      tr.probsg <- SampleDirichlet(tcountsg[gc, ] + 1)
      tr.probsg <- tr.probsg / sum(tr.probsg)
      sv <- which(lookup == gn)
      if (sum(sv)==0) {
        temp <- rep(1e-10, length(t2))
      } else {
        temp <- lpi[[sv]]
      }
      tr.probs[j-1] <- (tr.probsg[gn] * temp[ln]) / length(temp)
      temp.mat <- matrix(0, ncol=t2, nrow= t1)
      for(kk in 1:t2) {
        temp.mat[ , kk] <- rep((tr.probsg[gn] * temp[kk]) , t1)
      }
      temp.mat <- t(sweep(t(temp.mat), 2, rowSums(temp.mat), "/"))
      transMat[[j-1]] <- temp.mat
    }
    u[j-1] <- runif(1, 0, tr.probs[j-1])
  }
  ret <- list()
  ret$trans <- tr.probs
  ret$transMat <- transMat
  ret$u <- u
  return(ret)
}