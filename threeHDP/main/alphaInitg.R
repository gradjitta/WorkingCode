alphaInitg <- function(paths, L, lookup, phi, K0) {
  emit <- paths$emit
  Np <- length(paths$grids)
  N0 <- dim(phi)[2] / 2
  n1 <- which(paths$grids[[1]] == lookup)
  alpha <- lapply(1:Np, function(i) NULL)
  alpha[[1]] <- rep(1, K0)
  for (k in 1:K0) {
    if (sum(n1) == 0) {
      m <- 1
    } else {
      m <- k
    }
    alpha[[1]][k] <- prod(dnorm(emit[[1]], phi[m, 1:5], sqrt(phi[m, 6:10])))
  }
  for (i in 2:Np) {
    nc <- which(paths$grids[i] == lookup)
    if (sum(nc) == 0) {
      alpha[[i]] <- rep(0, 1)
    } else {
      alpha[[i]] <- rep(0, K0)
    }
  }
  return(alpha)
}
