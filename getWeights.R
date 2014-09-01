getWeights <- function(z) {
  rz <- c(1, cumprod(1 - z))[1:length(z)]
  w <- rz * z
  return (w)
}