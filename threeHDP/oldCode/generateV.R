generateV <- function(d, k, g, l, Ck) {
  v.temp <- rep(0, l)
  for (s in 1:l) {
    v.temp <- sampleVk(d, s+k, g, Ck)
  }
  return(v.temp)
}