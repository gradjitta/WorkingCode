updateL <- function(L, keep) {
  N <- length(L)
  for (i in 1:N) {
    temp <- rep(0,length(L[[i]]))
    for (j in 1:length(keep)) {
      temp <- temp + (L[[i]] == keep[j])*which(keep == keep[j])
    }
    L[[i]] <- temp
  }
  L
}
