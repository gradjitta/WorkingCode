getTransProbsfcp <- function(tcounts, cw, at, Kp ) {
  M <- 520
  pcounts1 <- tcounts[[1]]
  pcounts2 <- tcounts[[2]]
  pcountsg <- tcounts[[3]]
  time1 <- Sys.time()
  for (i in 1:M) {
    J <- length(tcounts[[1]][[i]])
    for (j in 1:J) {
      if (!is.null(tcounts[[1]][[i]][[j]])) {
        # K <-dim(tcounts[[1]][[i]][[j]])[1]
        tempMat1 <- (as.matrix(tcounts[[1]][[i]][[j]]))
        padding <- rep(0,Kp-ncol(tempMat1))
        gns <- getNbGrids(i)[getNbGrids(i) > 0]
        pcounts1[[i]][[j]] <- t((apply(tempMat1,1,function(x) rdirichlet (1,c(x,padding)+ cw[gns[j], ]*at  ))))
        if (dim(tempMat1)[2] == 1 && (Kp != 2 )) {
          pcounts1[[i]][[j]] <- t(pcounts1[[i]][[j]])
        }
      }
    }
    K <-dim(tcounts[[2]][[i]])[1]
    if (!is.null(K)) {
      tempMat2 <- (as.matrix(tcounts[[2]][[i]]) + 1)
      pcounts2[[i]] <- t((apply(tempMat2,1,function(x) rdirichlet(1,x))))
      if (dim(tempMat2)[2] == 1) {
        pcounts2[[i]] <- t(pcounts2[[i]])
      }
    }
  }
  tempMat3 <- (as.matrix(tcounts[[3]]) + 1)
  pcountsg <- t((apply(tempMat3,1,function(x) colMeans(rdirichlet(5,x)))))
  ret <- list()
  ret <- list(pcounts1, pcounts2,pcountsg)
  print(Sys.time() - time1)
  return(ret)
}

