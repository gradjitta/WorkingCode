getTransgnp <- function(g, pcounts1, pcounts2, pcountsg, L, K0, lpi, lookup, Kp, at, cSize) {
  #getTrans(paths[[i]]$grids, tcounts1, tcounts2, tcountsg,
  #                             L[[i]], sample.hdp$CK, sample.hdp$lpi, lookup)
  N <- length(g)
  MC <- 1e-100
  Lmin <- 100
  Bmin <- 80
  #tempMat3 <- (as.matrix(pcountsg) + 1)
  #pcountsg <- t((apply(tempMat3,1,function(x) colMeans(rdirichlet(10,x)))))
  tr.probs <- rep(0,(N-1))
  transMat <- lapply(1:(N-1), function(i) NULL)
  u <- rep(0, N-1)
  for(j in 2:N) {
    gn <- g[j]
    ln <- L[j]
    lc <- L[j-1]
    gc <- g[j-1]
    if( any(getNbGridsF(gc, cSize, Lmin, Bmin) == gn) ) {
      gk <- which(getNbGridsF(gc, cSize, Lmin, Bmin) == gn)
      t1.pmat <- pcounts1[[gc]][[gk]]
      t2.pmat <- pcounts2[[gc]]
      gk <- which(getNbGridsF(gc, cSize, Lmin, Bmin) == gn)
      temp.mat <- t1.pmat * t2.pmat[, gk]
      #######################
      kdiff <- Kp - K0
      mat.temp <- matrix(0, kdiff, Kp)
      if (kdiff > 0) {
        for (gind in 1:kdiff) {
           padding <- rep(0,kdiff)
           #trans1 <- rdirichlet(1,c(rep(1,K0),padding) + at*lpi[gn, ])
           trans1 <- rdirichlet(1,c(rep(0,K0),padding) + at*lpi[gn, ])
           gnb <- getNbGridsF(gc, cSize, Lmin, Bmin)
           trans2 <- rdirichlet(1,c(rep(0,length(gnb))) + 1)
           mat.temp[gind, ] <- trans1 * trans2[gk]
        }
      }
      temp.mat <- rbind(temp.mat, mat.temp)
      #######################
      #temp.mat <- t(sweep(t(temp.mat), 2, rowSums(temp.mat), "/"))
      transMat[[j-1]] <- temp.mat / rowSums(temp.mat)
      #tr.probs[j-1] <- temp.mat[lc, ln]
      tr.probs[j-1] <- transMat[[j-1]][lc, ln]
      #print(Sys.time() - time2)
    } else {
      kdiff <- Kp - K0
      tr.probsg <- pcountsg[gc, ]
      sv <- which(lookup == gn)
      if (sum(sv)==0) {
        #temp <- rep(1e-10, Kp)
        temp <- lpi[gn, ]
      } else {
        #temp <- lpi[sv, ]
         temp <- lpi[gn, ]
      }
      ###
      #temp <- rep(1/Kp,Kp)
      ###
      # tr.probs[j-1] <- (tr.probsg[gn] * temp[ln]) / length(temp)
      tr.probs[j-1] <- (tr.probsg[gn] * temp[ln])
      # - - - - - - - - - - - - - -
      temp.mat <- matrix(0, ncol=Kp, nrow= Kp)
      for(kk in 1:Kp) {
        temp.mat[ , kk] <- rep((tr.probsg[gn] * temp[kk]) , Kp)
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
      tr.probs[j-1] <- transMat[[j-1]][lc, ln]
      if (identical(transMat[[j-1]][lc, ln], numeric(0))) {
        print(paste(ln,lc,j))
      }
    }
    u[j-1] <- runif(1, 0, tr.probs[j-1])
    #print(paste(t1, t2))
  }
  ret <- list()
  ret$trans <- tr.probs
  ret$transMat <- transMat
  ret$u <- u
  return(ret)
}

