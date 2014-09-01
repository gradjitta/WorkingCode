hdpHmm <- function(paths, numIter, L, hdp.data, sample.hdp, Doom) {
  # paths1 <-  readRDS("Desktop/paths.rds")
  # saveRDS(paths, "Desktop/paths.rds")
  library(Matrix)
  # HDP HMM, using Beam Sampling to restrict the transitions and 
  # forward filtering backward sampling(FFBS) to train the HMM
  # 
  # Args:
  #   paths  : all the paths for the month of october used.
  #            Further each element in path has the following data
  #             -- samples or feature vectors corresponding to the path
  #             -- grid blocks visited
  #             -- states which is randomly assigned given the initial
  #                hdp states
  # 
  #   numIter : number of iterations
  # 
  # Returns:
  #  outputs from FFBS
  #    -- transition probabilities, local and global clusters/states 
  ################################################
  # 1) Initialize the transition matrix, 
  #    this is the ouput from a long HDP run lhssmmNEF.R with seed 4673
  #sample.hdp <- lhssmmNEF(1000, 0, 0, hdp.data)
  d <- sample.hdp$d
  bk <- sample.hdp$beta.k
  states <- sample.hdp$CK
  phi <- sample.hdp$phi
  # Extract paths, consisting of grid block and corresponding 
  # feature vector and constrained transition variable.
  # paths <- getPaths(hmmdata, sample.hdp$CK)
  iter <- 0
  M <-520
  G <- length(bk)
  # Initialize transition data structure
  lookup <- sapply(hdp.data[[2]], function(i) i)
  states <- sample.run1000$CK
  states <- Doom$Ck
  lookup <- Doom$lookup
  tcounts1 <- list() # P(L(t) | L(t-1), S(t))
  tcounts2 <- list() # P(S(t) | S(t-1), L(t-1))
  tcountsg <- Matrix(0, ncol= M, nrow = M) # for Jumps/ Missing transitions
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
          tcounts1[[i]][[j]] <- Matrix(0, ncol = length(states[[tj]]),
                                       nrow = length(states[[ti]]))
        } else {
          if(sum(ti) == 0 && sum(tj) != 0) {
              tcounts1[[i]][[j]] <- Matrix(0, ncol = length(states[[tj]]), nrow = 1)
            } else if(sum(ti) != 0 && sum(tj) == 0){
              tcounts1[[i]][[j]] <- Matrix(0, ncol = 1, nrow = length(states[[ti]]))
            } else {
              tcounts1[[i]][[j]] <- Matrix(0, ncol = 1, nrow = 1)
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
                                , nrow = length(states[[ti]]))
      } else {
        tcounts2[[i]] <- Matrix(0, ncol = length(near.i), nrow = 1)
      }
    }
  }
  # Init states
  # initL <- initStates(hdp.data[[2]], paths, sample.run1000$CK)
  # Main Loop
  iter <- 0
  
  while (iter < numIter) {
    # 1) reset counts
    if (iter > 0) {
      tcounts <- resetCounts(lookup, hdp.data[[2]], sample.run1000)
      tcounts1 <- tcounts[[1]]
      tcounts2 <- tcounts[[2]]
      tcountsg <- tcounts[[3]]
    }
    # Updates counts from paths
    bookL <- lapply(1:M, function(i) lapply(1:length(getNbGrids(i)), function(j) NULL ))
    for (i in 1:length(paths)) {
      print(i)
      if(length(paths[[i]]$grids) > 5) {
        for (j in 2:length(paths[[i]]$grids)) {
          gn <- paths[[i]]$grids[j]
          ln <- L[[i]][j]
          lc <- L[[i]][j-1]
          gc <- paths[[i]]$grids[j-1]
          gi <- which(lookup == gc)
          gk <- which(getNbGrids(gc) == gn)
          if( any(getNbGrids(gc) == gn) ) {
            tcounts1[[gc]][[gk]][lc, ln] <- tcounts1[[gc]][[gk]][lc, ln] + 1
            if ( class(bookL[[gc]][[gk]]) != "list") {
              if (sum(gi) != 0) {
                if ( (length(states[[gi]]) == 1) && (states[[gi]] == 1) ){ 
                  bookL[[gc]][[gk]] <- lapply(1:length(states[[gi]]), function(i) 1)
                } else {
                  bookL[[gc]][[gk]] <- lapply(1:length(states[[gi]]), function(i) NULL)
                }
              } else {
                bookL[[gc]][[gk]] <- lapply(1:1, function(i) NULL)
              }
            } else {
              bookL[[gc]][[gk]][[lc]] <- c(bookL[[gc]][[gk]][[lc]], ln)
            }
            tcounts2[[gc]][lc, gk] <- tcounts2[[gc]][lc, gk] + 1
          } else {
            tcountsg[gc, gn] <- tcountsg[gc, gn] + 1
          }
        }
      }
    }
    # Forward Filtering and Backward Sampling
    for (i in 1:length(paths)) {
      if(length(paths[[i]]$grids) > 5) {
        # 2) Calculate emission and transition probabilities
        trans.list <- getTrans(paths[[i]]$grids, tcounts1, tcounts2, tcountsg,
                               L[[i]], sample.hdp$CK, sample.hdp$lpi, lookup)
        trans.path <- trans.list$transMat
        emit.path <- getEmit(paths[[i]]$emit, paths[[i]]$grids, L[[i]],
                             sample.hdp$CK, sample.hdp$phi, lookup)
        # 3) get slices
        #runif(1, 0, getTrans(paths[[i]], j))
        u <- trans.list$u
        # 4) Forward Filtering
        # initialize alpha 
        alphas <- alphaInit(paths[[i]], L[[i]], sample.hdp$CK, lookup,
                            sample.hdp$phi) # returns a list of D(length of a path) elements
        D <- length(paths[[i]]$grids)
        for (t in 2:D) {
          temp.trans <- trans.path[[t-1]]
          sum.trans <- colSums(sweep((u[t-1] < temp.trans),MARGIN=1,alphas[[t-1]],`*`))
          if(emit.path[[t]] == 0) {
            temp.a <- (sum.trans)
          } else {
            temp.a <- (emit.path[[t]] * sum.trans)
          }
          if (sum(temp.a == 0)) {
            temp.a <- rep(1, length(sum.trans))
          }
          alphas[[t]] <- temp.a / sum(temp.a)
        }
        # 5) Backward Sampling
        # Initial backward sample
        #betas <-  betaInit(paths[[i]])
        # Sample the Dth stateQ
        L[[i]][D] <- 1 + sum(runif(length(alphas[[D]]),0,1) > cumsum(alphas[[D]]) )
        # Backward Sampling
        for (t in D:2) {
          temp.beta <-  alphas[[t-1]] * (u[t-1] < alphas[[t-1]])
          norm.beta <- temp.beta / sum(temp.beta)
          L[[i]][t-1] <- 1 + sum(runif(length(norm.beta),0,1) > cumsum(norm.beta) )
        }
      }
    }
    # 6) Sample phi, local sticks similar to HDP
    # L (but d for sample phi) from backward sampling used for the HDP updates or D
    dhdp <- dataForHdp(paths, L)
    data.hdp <- dhdp$data
    Lhdp <- dhdp$Lhdp
    d <- dhdp$d
    updates <- updateHDP(data, Ck, Lhdp, d, transition)
    bk <- updates$beta.k
    paths$states <- updates$CK
    phi <- updates$phi
    lpi <- updates$lpi
    iter <- iter + 1
  }
}