resetCounts <- function(lookup, hdp.data, sample.hdp){
  lookup <- sapply(hdp.data, function(i) i)
  states <- sample.hdp$CK
  tcounts1 <- list() # P(L(t) | L(t-1), S(t))
  tcounts2 <- list() # P(S(t) | S(t-1), L(t-1))
  tcountsg <- Matrix(0, ncol= M, nrow = M) # for Jumps/ Missing transitions
  for (i in 1:M) {
    if (any(c(1,26,495,520) == 1)) {
      near.i <- getNbGrids(i)
      for(j in 1:length(near.i)) {
        ti <- which(lookup == i)
        tj <- which(lookup == near.i[j])
        if(ti != 0 && tj != 0) {
          tcounts1[[i]][[j]] <- Matrix(0, ncol = length(states[[ti]]),
                                       nrow = length(states[[tj]]))
        } else {
          tcounts1[[i]][[j]] <- Matrix(0, ncol = 1, nrow = 1)
        }
      }
    }
  }
  for ( i in 1:M) {
    if (any(c(1,26,495,520) == 1)) {
      near.i <- getNbGrids(i)
      ti <- which(lookup == i)
      if(ti != 0) {
        tcounts2[[i]] <- Matrix(0, ncol = length(near.i) 
                                , nrow = length(states[[ti]]))
      } else {
        tcounts2[[i]] <- Matrix(0, ncol = length(near.i), nrow = 1)
      }
    }
  }
  ret <- list(tcounts1, tcounts2, tcountsg)
  return(ret)
}