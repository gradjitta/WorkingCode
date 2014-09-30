sampleVk <- function (d, k, g, Ck) {
  Ckp <- lapply(d, unique)
  nk <- sapply(Ckp, function(d) sum(d==k))
  a <- 1 + sum(nk)
  nnk <- sapply(Ckp, function(d) sum(d > k))
  b <- g + sum(nnk)
  z <- rbeta(1, a, b)
  return(z)
}