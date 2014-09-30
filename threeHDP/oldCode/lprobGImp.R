lprobGImp <- function(data, phi, L, j, m, J, i) {
  N <- dim(phi)[2] / 2
  ymap <- which(L[[j]] == i)
  y <- data[[j]][ymap, , drop =F]
  temp.m <- length(ymap) #sum(L[[j]] == i)
  #pd <- matrix(0, nrow = temp.m, ncol = N)
  pd <- rep(0,temp.m)
  for (feat in 1:N) {
    #pd[, feat] <- dnorm(y[, feat], phi[m, feat], sqrt(phi[m,N+feat]), log=TRUE)
    pd <- pd + dnorm(y[, feat], phi[m, feat], sqrt(phi[m,N+feat]), log=TRUE)
  }
  #pdl <- apply(pd,1,sum)
  #pdl <- rowSums(pd)
  if (length(pd) == 0) {
    pd <- -10000
  }
  return(pd)
}