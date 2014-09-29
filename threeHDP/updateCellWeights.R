updateCellWeights <- function(n,cellWeights,globalWeights,transAlpha,alpha,sN) {
  # Sample auxiliary variables telling to how many parts
  # each of the weights was split at the lower level
  m <- sampleTableCount(n,cellWeights,transAlpha,sN)

  # We only need their sums, padded with zeroes
  padding <- rep(1 ,length(globalWeights)-length(cellWeights))
  m <- c(colSums(m),padding)

  # Sample the actual weights given m, padding with zeroes
  cellWeights <- rdirichlet(1,m + alpha*globalWeights)

  return(list(cellWeights=cellWeights,m=m))
}
