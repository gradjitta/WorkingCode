lhssmmNEF <- function(num.iter, burn, trd, data) {
  set.seed(4673)
  J <- dim(summary(data))[1]
  initJ <- lapply(1:J,function(i) i)
  N0 <- dim(data[[1]])[2]
  d <- list()
  lc <- list(c(1,2,3), c(4,5,1), c(6, 7, 2))
  lc1 <- c(1,2,3, 4,5,1, 6, 7, 2)
  d <- lapply(initJ, function(i) d[[i]] <- rep(i, as.integer(dim(data[[i]])[1])))
  #d[[1]] <- c(rep(1, 50) , rep(2, 50), rep(3, 50))
  #d[[2]] <- c(rep(4, 50) , rep(5, 50), rep(1, 50))
  #d[[3]] <- c(rep(6, 50) , rep(7, 50), rep(2, 50))
  L <- list()
  L <- lapply(initJ, function(i) rep(1,dim(data[[i]])[1]))
  K <- lapply(d, unique)
  nks <- sapply(lapply(d,unique), max)
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
  lc <- list(c(1,2,3), c(4,5,1), c(6, 7, 2))
  #Ck <- lapply(lc, function(l) l)
  Ck <- lapply(initJ, function(l) 1)
  #Ck <- list(c(1,2,3))
  nkl <- sapply(lapply(Ck,unique), length)
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
  for (iter in 1:num.iter) {
    print("Iteration")
    print(iter)
    # Step 1, sample phi
    nkl <- sapply(lapply(Ck,unique), length)
    nks <- sapply(lapply(d, unique), max)
    njk <- max(nks)
    print("sampling phi")
    for (j in 1:njk) {
      ###################
      # sample mean with prior Normal(0, 1/s) on mean
      phi[j, 1:N0] <- sampleMeanjnm(data, d, j,
                                    phi[j,(N0+1):(2*N0),drop = T], s ,mean.init, J)
      # sample 1/variance with prior Gamma(epi, epi) on 1/variance
      phi[j, (N0+1):(2*N0)] <- sampleVariancejn(data,
                                                j, 1, d, phi[j, 1:N0, drop = T], J)
    }
    # Step 2, sample z(j)
    for (j in 1:njk) {
      v[j] <- sampleVk(d, j, gma, Ck)
    }
    beta.k <- getWeights(v)
    #########################
    # Step 3, sample global
    #
    # global slices
    for (j in 1:J) {
      for (l in 1:length(Ck[[j]])) {
        w[[j]][l] <- runif(1, 0, beta.k[Ck[[j]][l]])
      }
    }
    # Step 4, sample cjl (global assignments)
    print("generate new global sticks if needed")
    Ng <- 1 #length(beta.k)
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
    ## Change 1 start########################
    print("Calculate gprobs")
    for (j in 1:J) {
      gprobs[[j]] <- matrix(-10000, nrow = length(a[[j]]), ncol = length(beta.k))
      for (m in 1:length(beta.k)) {
        useLocalCluster <- which(beta.k[m] > w[[j]])
        for (i in useLocalCluster) {
          gprobs[[j]][i, m] <- sum(lprobGImp(data, phi, L, j, m, J, i))
        }
      }
    }
    ### Change 1 End#################################
    #print(gprobs)
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
    ###########################################
    # Local sticks and weights
    # Step 5, sample ajk
    for (j in 1:J) {
      for (l in 1:nkl[j]) {
        a[[j]][l] <- sampleAjk1(L, beta.k, alp, j, l, Ck)
      }
      lpi[[j]] <- getWeights(a[[j]])
    }
    # Step 6, sample local slices from uniform distribution
    for (j in 1:J) {
      for (l0 in 1:length(d[[j]])) {
        u[[j]][l0] <- runif(1, 0, lpi[[j]][L[[j]][l0]])
      }
    }
    # Step 7, sample djl ( local labels) from P(djl = l) where l is the local slice used
    print("Generate new Ck and local sticks if needed...")
    Nf <- rep(1,J)
    for (j in 1:J) {
      #print(j)
      while(TRUE) {
        if (sum(lpi[[j]][1:Nf[j]]) >= (1-min(u[[j]])) ) {
          break
        } else {
          Nf[j] <- Nf[j] + 1
          l <- Nf[j] - length(lpi[[j]])
          if (l > 0) {
            samplebk <- sample(length(beta.k), 1, prob=beta.k)
            Ck[[j]] <- c(Ck[[j]], samplebk)
            if (Nf[j] <= 49) {
              a[[j]] <- c(a[[j]], generateA(alp, l, d, beta.k, j, Ck, lpi, L))
            } else {
              a[[j]] <- c(a[[j]], 1)
              lpi[[j]] <- getWeights(a[[j]])
              break
            }
            lpi[[j]] <- getWeights(a[[j]])
          }
        }
      }
    }
#     print("Calc local probabilities")
#     for (j in 1:J) {
#       local.probs[[j]] <- matrix(0, nrow = dim(data[[j]])[1], ncol = length(lpi[[j]]) )
#       for (m in 1:length(lpi[[j]]) ) {
#         local.probs[[j]][, m] <- lprobT(data[[j]], lpi, phi, d, j, m, u, Ck)
#       }
#     }
#     for (l in 1:J) {
#       for (i in 1:dim(data[[l]])[1]) {
#         t.probs <- rep(0, length(local.probs[[l]][i, ]))
#         mapl <- which(local.probs[[l]][i, ] != 0)
#         t.probs[mapl] <- exp(local.probs[[l]][i, mapl])
#         L[[l]][i] <- sample(length(t.probs),1,prob=t.probs, replace=T)
#         d[[l]][i] <- Ck[[l]][L[[l]][i]]
#       }
#     }
    print("Calc local probabilities")
    for (j in 1:J) {
      local.probs[[j]] <- matrix(-10000, nrow = dim(data[[j]])[1], ncol = length(lpi[[j]]) )
      for (m in 1:length(lpi[[j]]) ) {
        local.probs[[j]][, m] <- lprobT(data[[j]], lpi, phi, d, j, m, u, Ck)
      }
    }
    for (l in 1:J) {
      #if(0) {
      for (i in 1:dim(data[[l]])[1]) {
        t.probs <- rep(-10000, length(local.probs[[l]][i, ]))
        t.probs <- exp(local.probs[[l]][i, ] - max(local.probs[[l]][i, ]))
        t.probs <- t.probs/sum(t.probs) #sum(exp(local.probs[[l]][i, ] - max(local.probs[[l]][i, ])))
        #mapl <- which(local.probs[[l]][i, ] != 0)
        #t.probs[mapl] <- exp(local.probs[[l]][i, mapl])
        L[[l]][i] <- sample(length(t.probs),1,prob=t.probs, replace=T)
        d[[l]][i] <- Ck[[l]][L[[l]][i]]
      }
      #}
      #rands <- runif(length(data[[l]])[1],0,1)
      #for (i in 1:dim(data[[l]])[1]) {
      #  t.probs <- rep(-10000, length(local.probs[[l]][i, ]))
      #  t.probs <- exp(local.probs[[l]][i, ] - max(local.probs[[l]][i, ]))
      #  t.probs <- t.probs/sum(t.probs) #sum(exp(local.probs[[l]][i, ] - max(local.probs[[l]][i, ])))
      #  t.probs <- cumsum(t.probs)
      #  L[[l]][i] <- length(which(rands[i]>t.probs)) + 1
      #  d[[l]][i] <- Ck[[l]][L[[l]][i]]
      #}
    }

    ####
    if(iter %% 51 == 0) {
      #print(d)
      #barplot(beta.k)
    }
    # Throw Ck
    print("throw Ck associated with empty clusters, sticks and weights associated with it")
    Ck <- mapply(function(i,j) Ck[[j]][1:i], lapply(L,max), initJ, SIMPLIFY = FALSE)
    a <- mapply(function(i,j) a[[j]][1:i], lapply(Ck,length), initJ, SIMPLIFY = FALSE)
    w <- mapply(function(i,j) w[[j]][1:i], lapply(Ck,length), initJ, SIMPLIFY = FALSE)
    for (j in 1:J) {
      lpi[[j]] <- getWeights(a[[j]])
    }
    # Throw phi
    print("Throw phi which was not used")
    dcut <-sapply(initJ, function(i) max(unique(Ck[[i]])) )
    phi <- phi[1:max(dcut), ,drop=FALSE]
    v <- v[1:max(dcut)]
    beta.k <- getWeights(v)
    print(d[[1]])
    print(paste("number of global sticks", length(beta.k)))
    #barplot(Nf)
    initL <- lapply(1:length(beta.k), function(i) i)
    somevals <- sapply(initL, function(i) sum(sapply(d, function(d) sum(d == i))))
    par(mfrow = c(1, 2), mar = c(1, 1, 1, 1), oma = c(0, 0, 2, 0), pty = "s")
    barplot(beta.k)
    barplot(somevals)
  }
  ret.val <- list()
  ret.val$CK <- Ck
  ret.val$d <- d
  ret.val$beta.k <- beta.k
  ret.val$lpi <- lpi
  ret.val$phi <- phi
  return(ret.val)
}