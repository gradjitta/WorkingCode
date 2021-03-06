pFFBSIhmmfnp <- function (i, paths, tcounts, initL, lookup, lpi, K0, phi, Kp, at) {
  set.seed(46)
  library(Matrix)
  tcounts1 <- tcounts[[1]]
  tcounts2 <- tcounts[[2]]
  tcountsg <- tcounts[[3]]
  L <- initL
  tempL <- L[[i]]
  D <- length(paths[[i]]$grids)
  if(D > 5) {
    # 2) Calculate emission and transition probabilities
    #time2 <- Sys.time()
    trans.list <- getTransgnp(paths[[i]]$grids, tcounts1, tcounts2, tcountsg,tempL, K0, lpi, lookup, Kp, at)
    trans.path <- trans.list$transMat
    emit.path <- getEmitg(Kp, paths[[i]]$emit, paths[[i]]$grids, phi, lookup)
    # 3) get slices
    #runif(1, 0, getTrans(paths[[i]], j))
    u <- trans.list$u
    # 4) Forward Filtering
    # initialize alpha 
    #time3 <- Sys.time()
    alphas <- alphaInitg(paths[[i]],tempL, lookup, phi, Kp)
    # returns a list of D(length of a path) elements
    #print(Sys.time() - time3)
    gpath <- paths[[i]]$grids
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
