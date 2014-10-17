alphaInitgc <- function(paths, L, lookup, phi, K0) {
  emit <- paths$emit
  Np <- length(paths$grids)
  N0 <- dim(phi)[2] / 2
  n1 <- which(paths$grids[[1]] == lookup)
  alpha <- lapply(1:Np, function(i) NULL)
  for (i in 1:Np) {
    nc <- which(paths$grids[i] == lookup)
    if (sum(nc) == 0) {
      alpha[[i]] <- rep(0, 1)
    } else {
      alpha[[i]] <- rep(0, K0)
    }
  }
  return(alpha)
}
