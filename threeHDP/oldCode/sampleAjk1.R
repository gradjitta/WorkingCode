sampleAjk1 <- function(L, beta.k, alpha, j, l, Ck) {
  a = 1 + sum(L[[j]]==l)
  b = 1 + sum(L[[j]] > l)
  ajk <- rbeta(1, a, b)
  return(ajk)
}