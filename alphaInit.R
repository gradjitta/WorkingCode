alphaInit <- function(paths, L, Ck, lookup, phi) {
  emit <- paths$emit
  Np <- length(paths$grids)
  N0 <- dim(phi)[2] / 2
  n1 <- which(paths$grids[[1]] == lookup)
  N <- length(Ck[[n1]])
  alpha <- lapply(1:Np, function(i) NULL)
  alpha[[1]] <- rep(1, N)
  for (k in 1:N) {
    if (sum(n1) == 0) {
      m <- 1
    } else {
      m <- Ck[[n1]][k]
    }
    pd <- 1
    for (feat in 1:N0) {
      pd <- pd * dnorm(emit[[1]][feat], phi[m, feat], sqrt(phi[m,N0+feat]))
    }
    alpha[[1]][k] <- pd
  }
  for (i in 2:Np) {
    nc <- which(paths$grids[[i]] == lookup)
    if (sum(nc) == 0) {
      alpha[[i]] <- rep(0, 1)
    } else {
      alpha[[i]] <- rep(0, length(Ck[[nc]]))
    }
  }
  return(alpha)
}