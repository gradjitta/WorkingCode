updateCountsNew <- function(M, gridlist, paths, lookup, states, L, tcounts) {
  corners <- c(1,26,495,520)
  #tcounts1 <- lapply(1:M, function(i) lapply(1:length(getNbGrids(i)), function(j) NULL ))
  #tcounts2 <- lapply(1:M, function(i) NULL)
  tcounts1 <- tcounts[[1]]
  tcounts2 <- tcounts[[2]]
  tcountsg <- matrix(0, ncol= M, nrow = M)
  time1 <- Sys.time()
  for (gcc in 1:M) {
    if (sum(corners == gcc) == 0) {
      temp.paths <- gridlist[[gcc]]
      for(i in temp.paths) {
        tl <- length(paths[[i]]$grids)
        gci <- which(paths[[i]]$grids == gcc)
        gll <- length(gci)
        gni <- gci + 1
        if(gni[gll] > tl) {
          gci <- gci[-gll]
          gni <- gni[-gll]
        }
        gns <- paths[[i]]$grids[gni]
        if (length(gns) != 0){
          for (j in 1:length(gns)) {
            gk <- which(getNbGrids(gcc) == gns[j])
            if( any(getNbGrids(gcc) == gns[j]) ) {
              if (!is.null(tcounts1[[gcc]][[gk]])) {
              ln <- L[[i]][gni[j]]
              lc <- L[[i]][gni[j]-1]
              tcounts1[[gcc]][[gk]][lc, ln] <- tcounts1[[gcc]][[gk]][lc, ln] + 1
              tcounts2[[gcc]][lc, gk] <- tcounts2[[gcc]][lc, gk] + 1
            }
            } else {
              tcountsg[gcc, gns[j]] <- tcountsg[gcc, gns[j]] + 1
            }
          }
        }
      }
      if (!is.null(gridlist[[gcc]])) {
        temp.iter <- length(getNbGrids(gcc))
        for (k in 1:temp.iter) {
          if (!is.null(tcounts1[[gcc]][[k]])) {
            tcounts1[[gcc]][[k]] <- Matrix(tcounts1[[gcc]][[k]], sparse=T)
          }
        }
        tcounts2[[gcc]] <- Matrix(tcounts2[[gcc]], sparse=T)
      }
    }
  }
  tcountsg <- Matrix(tcountsg, sparse=T)
  print(Sys.time() - time1)
  tcounts <- list(tcounts1, tcounts2, tcountsg)
  return(tcounts)
}
