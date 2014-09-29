updateGlobalWeights <- function(m,globalWeights,alpha,gamma,sN,maxNew=1) {
  # Sample auxiliary variables telling how often these were
  # split at the cell level
  l <- sampleTableCount(m,globalWeights,alpha,sN)

  # We only need their sums
  l <- colSums(l)

  # Sample the actual weights given m
  globalWeights <- rdirichlet(1,c(l,gamma))

  # Break the last weight into new sticks
  if(maxNew>1) {
    globalWeights <- globalWeights[1:(length(globalWeights)-1)]
    remaining <- 1-sum(globalWeights)
    betas <- rbeta(maxNew,1,gamma)
    # Explicit truncation after the last new one
    betas[maxNew] <- 1
    for(k in 1:maxNew) {
      globalWeights <- c(globalWeights,remaining*betas[k])
      remaining <- remaining*(1-betas[k])
    }
  }

  return(globalWeights)
}
