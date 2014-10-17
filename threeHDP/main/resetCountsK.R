resetCountsK <- function(M, K, lookup) {
  tcounts1 <- list() # P(L(t) | L(t-1), S(t))
  tcounts2 <- list() # P(S(t) | S(t-1), L(t-1))
  tcountsg <- matrix(0, ncol= M, nrow = M) # for Jumps/ Missing transitions
  corners <- c(1,26,495,520)
  tcounts1 <- lapply(1:M, function(i) lapply(1:length(getNbGrids(i)), function(j) NULL ))
  tcounts2 <- lapply(1:M, function(i) NULL)
  for (i in 1:M) {
    if (sum(corners == i) == 0) {
      near.i <- getNbGrids(i)
      for(j in 1:length(near.i)) {
        ti <- which(lookup == i)
        tj <- which(lookup == near.i[j])
        if(sum(ti) != 0 && sum(tj) != 0) {
          tcounts1[[i]][[j]] <- Matrix(0, ncol = K ,nrow = K)
        } else {
          if(sum(ti) == 0 && sum(tj) != 0) {
            tcounts1[[i]][[j]] <- Matrix(0, ncol = K, nrow = 1)
          } else if(sum(ti) != 0 && sum(tj) == 0){
            tcounts1[[i]][[j]] <- Matrix(0, ncol = K, nrow = K)
          } else {
            tcounts1[[i]][[j]] <- Matrix(0, ncol = K, nrow = 1)
          }
        }
      }
    }
  } 
  for ( i in 1:M) {
    if (any(c(1,26,495,520) == 1)) {
      near.i <- getNbGrids(i)
      
      ti <- which(lookup == i)
      if(sum(ti) != 0) {
        tcounts2[[i]] <- Matrix(0, ncol = length(near.i) 
                                , nrow = K)
      } else {
        tcounts2[[i]] <- Matrix(0, ncol = length(near.i), nrow = K)
      }
    }
  }
  retval <- list()
  retval <- list(tcounts1, tcounts2, tcountsg)
  return(retval)
}
