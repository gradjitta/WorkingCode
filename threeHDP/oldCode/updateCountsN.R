updateCountsN <- function(M, gridlist, paths, lookup, states, L) {
  corners <- c(1,26,495,520)
  tcounts1 <- lapply(1:M, function(i) lapply(1:length(getNbGrids(i)), function(j) NULL ))
  tcounts2 <- lapply(1:M, function(i) NULL)
  tcountsg <- matrix(0, ncol= M, nrow = M)
  time1 <- Sys.time()
  for (gc in 1:M) {
    if (sum(corners == gc) == 0) {
      temp.paths <- gridlist[[gc]]
      for(i in temp.paths) {
        tl <- length(paths[[i]]$grids)
        gci <- which(paths[[i]]$grids == gc)
        gll <- length(gci)
        gni <- gci + 1
        if(gni[gll] > tl) {
          gci <- gci[-gll]
          gni <- gni[-gll]
        }
        gns <- paths[[i]]$grids[gni]
        if (length(gns) != 0){
        for (j in 1:length(gns)) {
          gk <- which(getNbGrids(gc) == gns[j])
          if( any(getNbGrids(gc) == gns[j]) ) {
            if (is.null(tcounts1[[gc]][[gk]])) {
              ti <- which(lookup == gc)
              tj <- which(lookup == gns[j])
              tcounts1[[gc]][[gk]] <- matrix(0, ncol = length(states[[tj]]), nrow = length(states[[ti]]))
              tcounts2[[gc]] <- matrix(0, ncol = length(getNbGrids(gc)), nrow = length(states[[ti]]))
              ln <- L[[i]][gni[j]]
              lc <- L[[i]][gni[j]-1]
              tcounts1[[gc]][[gk]][lc, ln] <- tcounts1[[gc]][[gk]][lc, ln] + 1
              tcounts2[[gc]][lc, gk] <- tcounts2[[gc]][lc, gk] + 1
              } else {
                ln <- L[[i]][gni[j]]
                lc <- L[[i]][gni[j]-1]
                tcounts1[[gc]][[gk]][lc, ln] <- tcounts1[[gc]][[gk]][lc, ln] + 1
                tcounts2[[gc]][lc, gk] <- tcounts2[[gc]][lc, gk] + 1
              }
          } else {
            tcountsg[gc, gns[j]] <- tcountsg[gc, gns[j]] + 1
          }
        }
        }
      }
      if (!is.null(gridlist[[gc]])) {
        temp.iter <- length(getNbGrids(gc))
        for (k in 1:temp.iter) {
          if (!is.null(tcounts1[[gc]][[k]])) {
            tcounts1[[gc]][[k]] <- Matrix(tcounts1[[gc]][[k]], sparse=T)
          }
        }
        tcounts2[[gc]] <- Matrix(tcounts2[[gc]], sparse=T)
      }
    }
  }
  tcountsg <- Matrix(tcountsg, sparse=T)
  print(Sys.time() - time1)
  tcounts <- list(tcounts1, tcounts2, tcountsg)
  return(tcounts)
}
