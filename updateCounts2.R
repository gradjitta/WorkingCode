updateCounts2 <- function(ii, M, paths, lookup, states, initcounts,L, spaths1,spaths2) {
  tcounts1 <- initcounts[[1]]
  tcounts2 <- initcounts[[2]]
  tcountsg <- initcounts[[3]]
  bookL <- lapply(1:M, function(i) lapply(1:length(getNbGrids(i)), 
                                           function(j) if (any(lookup == i)) {lapply(1:length(states[[which(lookup == i)]]), 
                                                                                     function(p) NULL)} else {NULL}))
  # bookL <- lapply(1:M, function(i) lapply(1:length(getNbGrids(i)), function(j) NULL ))
  for (i in spaths1[ii]:spaths2[ii]) {
    if(length(paths[[i]]$grids) > 5) {
      for (j in 2:length(paths[[i]]$grids)) {
        gn <- paths[[i]]$grids[j]
        ln <- L[[i]][j]
        lc <- L[[i]][j-1]
        gc <- paths[[i]]$grids[j-1]
        gi <- which(lookup == gc)
        gk <- which(getNbGrids(gc) == gn)
        if( any(getNbGrids(gc) == gn) ) {
          tcounts1[[gc]][[gk]][lc, ln] <- tcounts1[[gc]][[gk]][lc, ln] + 1
          #bookL[[gc]][[gk]][[lc]][length(which(bookL[[gc]][[gk]][[lc]] > 0)) + 1] = ln
          #bookL[[gc]][[gk]][[lc]] <- c(bookL[[gc]][[gk]][[lc]], ln)
          tcounts2[[gc]][lc, gk] <- tcounts2[[gc]][lc, gk] + 1
        } else {
          tcountsg[gc, gn] <- tcountsg[gc, gn] + 1
        }
      }
    }
  }
 return(list(tcounts1,tcounts2,tcountsg,bookL))
}
