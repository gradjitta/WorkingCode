updateLpic <- function(L, bookL, beta.k, Ck, data, tcounts, lookup) {
  J <- length(lookup)
  #lapply( 1:J, function(i) lapply(1:length(getNbGrids(lookup[i])), function(j) if (sum(which(getNbGrids(lookup[i])[j] == lookup)) > 0)
  #                                 lapply( 1:length(Ck[[getNbGrids(lookup[i])[j] == lookup]]), function(k) NULL )))
  #
  corners <- c(1,26,495,520)
  tempCk <- Ck
  Nf <- rep(1, J)
  newNF <- rep(1, J)
  lpi <- lapply(1:J, function(i) lapply(1:length(getNbGrids(lookup[i])), function(j)  lapply(1:length(Ck[[i]]), function(k) NULL) ))
  a <- lapply(1:J, function(i) lapply(1:length(getNbGrids(lookup[i])), function(j)  lapply(1:length(Ck[[i]]), function(k) NULL) ))
  u <- lapply(1:J, function(i) lapply(1:length(getNbGrids(lookup[i])), function(j)  lapply(1:length(Ck[[i]]), function(k) NULL) ))
  for (j in lookup) {
    nbg <- getNbGrids(j)
    for (i in 1:length(nbg)) {
      if ((sum(corners == i) == 0) && (sum(nbg[i] == lookup) != 0) && (!is.null(bookL[[j]][[i]])) ) {
        ff <- which(lookup == j)
        fn <- which(lookup == nbg[i])
      tl <- tempCk[[ff]]
      for (k in 1:length(tl)) {
        temp.l <- tcounts[[j]][[i]][k, ] + 1
        for(l in 1:length(temp.l)){
          a[[ff]][[i]][[k]] <- c(a[[ff]][[i]][[k]], isampleAjk1(temp.l, beta.k, alp, j, l, Ck))
        }
        lpi[[ff]][[i]][[k]] <-  getWeights(a[[ff]][[i]][[k]])
        # sample local slices
        #tL <- bookL[[j]][[i]][[k]]
        #if (is.null(tL)) {
        #  tL <- 1
        #}
        # change
        tL <- rep(which(tcounts[[j]][[i]][k, ] > 0), tcounts[[j]][[i]][k, ][tcounts[[j]][[i]][k, ] > 0])
        if (sum(tL) == 0) {
           tL <- 1
        }
        for (l0 in 1:length(tL)) {
          u[[ff]][[i]][[k]][l0] <- runif(1, 0, lpi[[ff]][[i]][[k]][tL[l0]])
        }
        tlpi <- lpi[[ff]][[i]][[k]]
        tu <- u[[ff]][[i]][[k]]
        ta <- a[[ff]][[i]][[k]]
        Nfj <- length(tlpi)
        # Generate new Ck and local sticks if needed...
        while(TRUE) {
          if (sum(tlpi[1:Nfj]) >= (1-min(tu)) ) {
            break
          } else {
            newNfj <- length(Ck[[fn]])
            Nfj <- Nfj + 1
            l <- Nfj - length(tlpi)
            l1 <- Nfj -newNfj
            if (l > 0) {
              if (l1 > 0) {
                samplebk <- sample(length(beta.k), 1, prob=beta.k)
                Ck[[fn]] <- c(Ck[[fn]], samplebk)
              }
              if (newNfj <= 49) {
                ta <- c(ta, isampleAjk1(temp.l, beta.k, 1, j, (1+length(tlpi)), Ck))
              } else {
                ta <- c(ta, 1)
                tlpi <- getWeights(ta)
                lpi[[ff]][[i]][[k]] <- tlpi
                a[[ff]][[i]][[k]] <- ta
                break
              }
              tlpi <- getWeights(ta)
              lpi[[ff]][[i]][[k]] <- tlpi
              a[[ff]][[i]][[k]] <- ta
            }
          }
        }
      }
    }
  }
  }
  local.params <- list()
  local.params$Ck <- Ck
  local.params$lpi <- lpi
  local.params$a <- a
  return(local.params)
}

