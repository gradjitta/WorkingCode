initStates <- function(lookup, paths, Ck) {
  #lookup <- sapply(hdpc, function(i) i)
  L <- list()
  for(i in 1:length(paths)) {
    g <- paths[[i]]$grids
    tl <- rep(0, length(g))
    for (j in 1:length(g)) {
      tc <- which(lookup == g[j])
      if(length(tc) != 0) {
        vals <- runif(length(Ck[[tc]]))
        sampC <- which(vals == max(vals))
        tl[j]  <- sampC
      } else {
        tl[j]  <- 1
      }
    }
    L[[i]] <- tl
  }
  return(L)
}
