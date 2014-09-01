getTransProbs <- function(tcounts) {
  pcounts1 <- tcounts[[1]]
  pcounts2 <- tcounts[[2]]
  pcountsg <- tcounts[[3]]
  time1 <- Sys.time()
  for (i in 1:520) {
    J <- length(tcounts[[1]][[i]])
      for (j in 1:J) {
        if (!is.null(tcounts[[1]][[i]][[j]])) {
         K <-dim(tcounts[[1]][[i]][[j]])[1]
          for (k in 1:K) {
            #if (sum(tcounts[[1]][[i]][[j]][k, ]) != 0) {
            pcounts1[[i]][[j]][k, ] <- SampleDirichlet(tcounts[[1]][[i]][[j]][k, ] + 1)
          #}
          }
        }
      }
  }
  for (i in 1:520) {
    K <-dim(tcounts[[2]][[i]])[1]
    for (k in 1:K) {
      pcounts2[[i]][k, ] <- SampleDirichlet(tcounts[[2]][[i]][k, ] + 1)
    }
  }
  for (i in 1:520) {
    pcountsg[i, ] <- SampleDirichlet(tcounts[[3]][i , ] + 1)
  }
  print(Sys.time() - time1)
  ret <- list()
  ret <- list(pcounts1, pcounts2,pcountsg)
  return(ret)
}
