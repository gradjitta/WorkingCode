pFFBSIhmmf <- function (i, paths, tcounts, initL, sample.run1000, lookup) {
  set.seed(46)
  library(Matrix)
  tcounts1 <- tcounts[[1]]
  tcounts2 <- tcounts[[2]]
  tcountsg <- tcounts[[3]]
  L <- initL
  tempL <- L[[i]]
  sample.hdp <- sample.run1000
  D <- length(paths[[i]]$grids)
  if(D > 5) {
    # 2) Calculate emission and transition probabilities
    #time2 <- Sys.time()
    trans.list <- getTransp(paths[[i]]$grids, tcounts1, tcounts2, tcountsg,
                           tempL, sample.hdp$CK, sample.hdp$lpi, lookup)
    trans.path <- trans.list$transMat
    emit.path <- getEmitf(paths[[i]]$emit, paths[[i]]$grids, tempL,
                         sample.hdp$CK, sample.hdp$phi, lookup)
    # 3) get slices
    #runif(1, 0, getTrans(paths[[i]], j))
    u <- trans.list$u
    # 4) Forward Filtering
    # initialize alpha 
    #time3 <- Sys.time()
    alphas <- alphaInit(paths[[i]], tempL, sample.hdp$CK, lookup,
                        sample.hdp$phi) # returns a list of D(length of a path) elements
    #print(Sys.time() - time3)
    for (t in 2:D) {
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
    # 5) Backward Sampling
    # Initial backward sample
    #betas <-  betaInit(paths[[i]])
    # Sample the Dth stateQ
    #time4 <- Sys.time()
    tempL[D] <- 1 + sum(runif(length(alphas[[D]]),0,1) > cumsum(alphas[[D]]) )
    # Backward Sampling
    for (t in D:2) {
      temp.pi <- trans.path[[t-1]][ , tempL[t]]
      temp.beta <-  alphas[[t-1]] * (u[t-1] < temp.pi)
      if(sum(temp.beta) == 0) {
        temp.beta <- alphas[[t-1]]
      }
      norm.beta <- temp.beta / sum(temp.beta)
      tempL[t-1] <- 1 + sum(runif(length(norm.beta),0,1) > cumsum(norm.beta) )
    }
    if (sum(is.na(tempL)) > 0) {
      print("have NAs")
    }
    #print(Sys.time() - time2)
  }
  return(tempL)
}

