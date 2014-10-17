getMCounts <- function(M, tcounts1, keep) {
K <- length(keep)
m <- list()
for(i in 1:M) {
  m[[i]] <- matrix(rep(0,K), 1, K)
  if (any(c(1,26,495,520) == i)) {
    m[[i]] <- rbind(m[[i]],rep(0,K))
  } else {
    L1 <- sapply(getNbGrids(i),function(p)  which(getNbGrids(p) == i))
    L2 <- getNbGrids(i)
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
