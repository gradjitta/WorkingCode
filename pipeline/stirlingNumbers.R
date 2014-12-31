# Precompute (scaled) stirling numbers up to N,N
#

stirlingNumbers <- function(N) {
  res <- matrix(0,N+1,N+1)
  res[1,1] <- 1
  for(n in 1:N) {
    for(k in 1:n) {
      res[n+1,k+1] <- -(n-1)*res[n,k+1] + res[n,k]
    }
    # Scale
    res[n+1,] <- res[n+1,] / max(res[n+1,])
  }
  return(abs(res[2:(N+1),2:(N+1)]))
}
