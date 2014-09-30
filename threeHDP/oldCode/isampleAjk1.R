isampleAjk1 <- function(L, beta.k, alpha, j, l, Ck) {
  if(length(L) < l) {
    a = 1
    b = 1
  } else {
    a = 1 + L[l]
    b = 1 + sum(L[-l])
  }
  ajk <- rbeta(1, a, b)
  return(ajk)
}