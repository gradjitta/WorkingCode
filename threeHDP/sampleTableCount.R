sampleTableCount <- function(n,beta,alpha,sN) {
  m <- 0*n

  for(i in 1:nrow(n)) {
    for(k in 1:length(beta)) {
      if(n[i,k]>0) {
        prob <- sN[n[i,k],1:n[i,k]]*(alpha*beta[k])^(1:n[i,k])
        prob <- prob/sum(prob)
        m[i,k] <- sample(n[i,k],1,prob=prob)
      } else {
        m[i,k] <- 1
      }
    }
  }
  return(m)
}
