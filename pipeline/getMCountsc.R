getMCountsc <- function(M, tcounts1, keep,cSize) {
K <- length(keep)
Lmin <- 100
Bmin <- 80
corners <- getCorners(cSize, Lmin, Bmin)
m <- list()
keep <- 1:K
for(i in 1:M) {
  m[[i]] <- matrix(rep(0,K), 1, K)
  if (any(corners == i)) {
    m[[i]] <- rbind(m[[i]],rep(0,K))
  } else {
    L1 <- sapply(getNbGridsF(i, cSize, Lmin, Bmin),function(p)  which(getNbGridsF(p, cSize, Lmin, Bmin) == i))
    L2 <- getNbGridsF(i, cSize, Lmin, Bmin)
    tempm <- rep(0, K)
    for (ind in 1:length(L2)) {
        ii <- L2[ind]
        j <- L1[ind]
        if(sum(dim(tcounts1[[ii]][[j]])) == 2 && K != 1 ) {
          m[[i]] <- rbind(m[[i]],rep(0,K))
        } else {
          if(!is.null(tcounts1[[ii]][[j]])) {
            tempm <- as.matrix(tcounts1[[ii]][[j]])
            tempm <- tempm[, keep, drop = F]
            m[[i]] <- rbind(m[[i]], tempm)
          } else {
            m[[i]] <- rbind(m[[i]], rep(0,K))
          }
        }
    }
  }
  m[[i]] <- m[[i]][-1, , drop = F]
}
m
}

