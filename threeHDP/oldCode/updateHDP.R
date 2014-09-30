updateHDP <- function(data, Ck, L, d, lookup, tcounts, bookL) {
  # HDP step to update phi, Ck, local sticks and global sticks
  # 
  # Args:
  #   data: paths$samples
  #   Ck: local cluster to global cluster 
  #   L : the local cluster associated with each sample generated from
  #       backward sampling step in HDPHMM
  #   transition: the complete transition matrix which needs to be 
  #               updated based on the number of local and global clusters
  #
  # Returns:
  #  Updated phi, Ck, local states and transition matrices
  ################################################
  # Initialize 
  J <- dim(summary(data))[1]
  initJ <- lapply(1:J,function(i) i)
  N0 <- dim(data[[1]])[2]
  K <- lapply(d, unique)
  nks <- sapply(lapply(Ck,unique), max)
  njk <- max(nks)
  v <- rep(0, njk)
  gma <- 1 # gamma parameter for global stick breaking process
  alp <- 1 # alpha parameter for local ,,,,
  #init local sticks
  a <- lapply(1:J, function(i) 0)
  lpi <- lapply(1:J, function(i) 0)
  u <- lapply(1:J, function(i) 0)
  #init global stick weights
  w <- lapply(1:J, function(i) 0)
  #Ck <- lapply(lc, function(l) l)
  #Ck <- lapply(initJ, function(l) 1)
  #Ck <- list(c(1,2,3))
  nkl <- sapply(1:J, function(i) length(Ck[[i]]))
  gprobs <- lapply(1:J, function(i) matrix(0,1,1))
  local.probs <- lapply(1:J, function(i) matrix(0,1,1))
  mean.init <- 1.5
  s <- 1
  phi <- matrix(0, njk, 2*N0)
  for (i in 1:njk) {
    m0 <- rnorm(N0, mean.init, sqrt(1/s))
    v0 <- 1/rgamma(N0, 1, 1)
    phi[i, ] <- c(m0, v0)
  }
  for (j in 1:njk) {
    v[j] <- sampleVk(d, j, gma, Ck)
  }
  beta.k <- getWeights(v)
  for (j in 1:J) {
    for (l in 1:nkl[j]) {
      a[[j]][l] <- sampleAjk1(L, beta.k, alp, j, l, Ck)
    }
    lpi[[j]] <- getWeights(a[[j]])
  }
  nks <- sapply(lapply(Ck, unique), max)
  njk <- max(nks)
  #
  print("Sample phi..")
  for (j in 1:njk) {
    ###################
    # sample mean with prior Normal(0, 1/s) on mean
    phi[j, 1:N0] <- sampleMeanjnm(data, d, j,
                                  phi[j,(N0+1):(2*N0),drop = T], s ,mean.init, J)
    # sample 1/variance with prior Gamma(epi, epi) on 1/variance
    phi[j, (N0+1):(2*N0)] <- sampleVariancejn(data,
                                              j, 1, d, phi[j, 1:N0, drop = T], J)
  }
  # Sample global sticks
  for (j in 1:njk) {
    v[j] <- sampleVk(d, j, gma, Ck)
  }
  beta.k <- getWeights(v)
  #
  # global slices
  for (j in 1:J) {
    for (l in 1:length(Ck[[j]])) {
      w[[j]][l] <- runif(1, 0, beta.k[Ck[[j]][l]])
      if (is.na(w[[j]][l])) {
         print(beta.k)
         print(Ck[[j]])
         print(l)
      }
    }
  }
  # generate new global sticks if needed
  print("generate new global sticks if needed")
  Ng <- 1
  minw <- min(sapply(1:J, function(j) min(w[[j]]) ))
  while(TRUE) {
      if ((sum(beta.k[1:Ng])) >= (1-minw) ) {
        break
      } else {
        if (length(beta.k) == 400) {
          v <- c(v[1:399], 1)
          beta.k <- getWeights(v)
          break
        }
        Ng <- Ng + 1
        l <- Ng - length(beta.k)
        if(Ng == 400) {
          print(Ng)
          print(length(beta.k))
          if (l > 0) {
            v <- c(v, 1)
            beta.k <- getWeights(v)
            ms <- rnorm(N0, mean.init, sqrt(1/s))
            vs <- 1/rgamma(N0, 1, 1)
            new.mat <- matrix(c(ms,vs), 1, 2*N0)
            phi <- rbind(phi, new.mat)
          }
        } else {
          if (l > 0) {
            v <- c(v, generateV(d, length(beta.k), gma, l, Ck))
            beta.k <- getWeights(v)
            ms <- rnorm(N0, mean.init, sqrt(1/s))
            vs <- 1/rgamma(N0, 1, 1)
            new.mat <- matrix(c(ms,vs), 1, 2*N0)
            phi <- rbind(phi, new.mat)
          }
        }
      }
  }
  # Calculate gprobs
    for (j in 1:J) {
      gprobs[[j]] <- matrix(-10000, nrow = length(a[[j]]), ncol = length(beta.k))
      for (m in 1:length(beta.k)) {
        useLocalCluster <- which(beta.k[m] > w[[j]])
        for (i in useLocalCluster) {
          gprobs[[j]][i, m] <- sum(lprobGImp(data, phi, L, j, m, J, i))
        }
      }
    }
  # Updating Ck ...
  print("Updating Ck ...")
  for (j in 1:J) {
      for (i in 1:length(Ck[[j]])) {
        # Check!
        temp.probs <- rep(0, length(gprobs[[j]][i, ]))
        temp.probs <- exp(gprobs[[j]][i, ] - max(gprobs[[j]][i, ])) / sum(exp(gprobs[[j]][i, ] - max(gprobs[[j]][i, ])))
        Ck[[j]][i] <- sample(1:length(gprobs[[j]][i, ]),1 
                             ,prob=temp.probs, replace=T)
      }
  }
  # Local sticks and weights
  # sample local sticks(L, bookL, beta.k, Ck, data, tcounts, lookup)
  # (L, bookL, beta.k, Ck, data, tcounts, lookup)
  # updateLpi(L, bookL, new1000$beta.k, DataHdp$Ck, DataHdp$data, tcounts[[1]], DataHdp$lookup)
  local.params <- updateLpic(L, bookL, beta.k, Ck, data, tcounts[[1]], lookup)
  Ck <- local.params$Ck
  # Throw Ck
  print("throw Ck associated with empty clusters, sticks and weights associated with it")
  Ck <- mapply(function(i,j) Ck[[j]][1:i], lapply(L,max), initJ, SIMPLIFY = FALSE)
  a <- mapply(function(i,j) a[[j]][1:i], lapply(Ck,length), initJ, SIMPLIFY = FALSE)
  w <- mapply(function(i,j) w[[j]][1:i], lapply(Ck,length), initJ, SIMPLIFY = FALSE)
  for (j in 1:J) {
    lpi[[j]] <- getWeights(a[[j]])
  }
  # Throw phi
  dcut <-sapply(initJ, function(i) max(unique(Ck[[i]])) )
  phi <- phi[1:max(dcut), ,drop=FALSE]
  v <- v[1:max(dcut)]
  beta.k <- getWeights(v)
  # Update transition matrix
  ret.val <- list()
  ret.val$CK <- Ck
  ret.val$beta.k <- beta.k
  ret.val$phi <- phi
  ret.val$lpi <- lpi
  return(ret.val)
  }
