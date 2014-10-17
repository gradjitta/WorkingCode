updateCountsS1 <- function(M, gridlist, paths, lookup, K, L, tcounts) {
  corners <- c(1,26,495,520)
  tcounts1 <- tcounts[[1]]
  tcounts2 <- tcounts[[2]]
  tcountsg <- matrix(0, ncol= M, nrow = M)
  jumpCounts <- rep(0, K)
  time1 <- Sys.time()
  for (gcc in 1:520) {
    NBS <- getNbGrids(gcc)
    tempMat2 <- as.matrix(tcounts2[[gcc]])
    tcounts1[[gcc]] <- lapply(1:length(tcounts1[[gcc]]), function(i)  if(!is.null(tcounts1[[gcc]][[i]])){as.matrix(tcounts1[[gcc]][[i]])})
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
            gk <- which(NBS == gns[j])
            if( any(NBS == gns[j]) ) {
              if (!is.null(tcounts1[[gcc]][[gk]])) {
                ln <- L[[i]][gni[j]]
                lc <- L[[i]][gni[j]-1]
                tcounts1[[gcc]][[gk]][lc, ln] <- tcounts1[[gcc]][[gk]][lc, ln] + 1
                tempMat2[lc, gk] <- tempMat2[lc, gk] + 1
              }
            } else {
              jumpCounts[L[[i]][gni[j]]] <- jumpCounts[L[[i]][gni[j]]] + 1
              tcountsg[gcc, gns[j]] <- tcountsg[gcc, gns[j]] + 1
            }
          }
        }
      }
      if (!is.null(gridlist[[gcc]])) {
        temp.iter <- length(NBS)
        for (k in 1:temp.iter) {
          if (!is.null(tcounts1[[gcc]][[k]])) {
            tcounts1[[gcc]][[k]] <- Matrix(tcounts1[[gcc]][[k]], sparse=T)
          }
        }
        tcounts2[[gcc]] <- Matrix(tempMat2, sparse=T)
      }
    }
  }
  tcountsg <- Matrix(tcountsg, sparse=T)
  print(Sys.time() - time1)
  tcounts <- list(tcounts1, tcounts2, tcountsg, jumpCounts)
  return(tcounts)
}
