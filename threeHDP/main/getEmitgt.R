getEmitgt <- function(K, emit, go, phi, lookup, N1) {
  # getEmit(paths[[i]]$emit, paths[[i]]$grids, L[[i]],
  #                         sample.hdp$CK, sample.hdp$phi, lookup)
  N0 <- length(go)
  N <- dim(phi)[2] / 2
  emit.probs <- lapply(1:N0, function(i) NULL)
  m <- 1:K
  for (i in 1:N0) {
    temp.emit <- rep(0, length(m))
    for(j in 1:length(m)) {
      mj <- m[j]
      temp.emit[j] <- sum(dnorm(emit[[i]], phi[mj, 1:N1], sqrt(phi[mj, (N1+1):(2*N1)]), log = TRUE))
    }
    num.temp <- exp(temp.emit - max(temp.emit))
    emit.probs[[i]] <- num.temp / sum(num.temp)
  }
  return(emit.probs)
}
