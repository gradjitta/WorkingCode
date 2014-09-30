getEmitf <- function(emit, go, L, Ck, phi, lookup) {
  # getEmit(paths[[i]]$emit, paths[[i]]$grids, L[[i]],
  #                         sample.hdp$CK, sample.hdp$phi, lookup)
  N0 <- length(go)
  N <- dim(phi)[2] / 2
  emit.probs <- lapply(1:N0, function(i) NULL)
  for (i in 1:N0) {
    if (sum(which(lookup == go[i])) == 0) {
      m <- 1 # This will never happen!!
    } else {
      gi <- which(lookup == go[i])
      m <- Ck[[gi]]
    }
    temp.emit <- rep(0, length(m))
    for(j in 1:length(m)) {
      mj <- Ck[[gi]][j]
      #pd <- 0
      #for (feat in 1:N) {
      #  pd <- pd + dnorm(emit[[i]][feat], phi[mj, feat], sqrt(phi[mj,N+feat]), log = TRUE)
      #}
      #temp.emit[j] <- pd 
      temp.emit[j] <- sum(dnorm(emit[[i]], phi[mj, 1:5], sqrt(phi[mj, 6:10]), log = TRUE))
    }
    num.temp <- exp(temp.emit - max(temp.emit))
    emit.probs[[i]] <- num.temp / sum(num.temp)
  }
  return(emit.probs)
}
