dataForHdp <- function(paths, L, Ck, lookup) {
  nCk <- list()
  N <- length(paths)
  data <- list()
  Lhdp <- list()
  d <- list()
  dataF <- list()
  LhdpF <- list()
  dF <- list()
  for (ini in 1:520) {
    data[[ini]] <- matrix(1:5, 1, 5)
    if(sum(ini == lookup) == 0) {
      nCk[[ini]] <- 1
    } else {
      nl <- which(lookup == ini)
      nCk[[ini]] <- Ck[[nl]]
    }
  }
  Lhdp <- lapply(1:520, function(i) c())
  d <- lapply(1:520, function(i) c())
  for(i in 1:N) {
    if(length(paths[[i]]$grids) > 5) {
    for(j in 1:length(paths[[i]]$grids)) {
      nk <- which(paths[[i]]$grids[j] == lookup)[1]
      data[[paths[[i]]$grids[j]]] <- rbind(data[[paths[[i]]$grids[j]]],
                                           paths[[i]]$emit[[j]])
      Lhdp[[paths[[i]]$grids[j]]]  <- c(Lhdp[[paths[[i]]$grids[j]]],
                                        L[[i]][j])
      if (!is.na(nk)) {
        d[[paths[[i]]$grids[j]]]  <- c(d[[paths[[i]]$grids[j]]],
                                       Ck[[nk]][L[[i]][j]])
      } else {
        d[[paths[[i]]$grids[j]]]  <- c(d[[paths[[i]]$grids[j]]], 1)
      }
      
    } 
    }
  }
  count <- 1
  for (i in 1:520) {
    if (dim(data[[i]])[1] != 1) {
      dataF[[count]] <- data[[i]][-1, ]
      LhdpF[[count]] <- Lhdp[[i]]
      dF[[count]] <- d[[i]]
      Ck[[count]] <- nCk[[i]]
      lookup[[count]] <- i
      count <- count + 1
    }
  }
  ret <- list()
  ret$d <- dF
  ret$Ck <- Ck
  ret$Lhdp <- LhdpF
  ret$data <- dataF
  ret$lookup <- lookup
  return(ret)
}