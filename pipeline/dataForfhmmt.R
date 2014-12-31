dataForfhmmt <- function(paths, L, lookup, N0, M) {
  N <- length(paths)
  Lhdp <- list()
  dataF <- list()
  LhdpF <- list()
  Lhdp <- lapply(1:M, function(i) c())
  data <- lapply(1:M, function(ini) matrix(1, nrow=1, ncol = N0) )
  d <- lapply(1:M, function(i) c())
  for(i in 1:N) {
    if(length(paths[[i]]$grids) > 5) {
      for(j in 1:length(paths[[i]]$grids)) {
        nk <- which(paths[[i]]$grids[j] == lookup)[1]
        data[[paths[[i]]$grids[j]]] <- rbind(data[[paths[[i]]$grids[j]]],  paths[[i]]$emit[[j]])
        Lhdp[[paths[[i]]$grids[j]]]  <- c(Lhdp[[paths[[i]]$grids[j]]], L[[i]][j])
      } 
    }
  }
  count <- 1
  for (i in 1:M) {
    if (dim(data[[i]])[1] != 1) {
      dataF[[count]] <- data[[i]][-1, , drop = F]
      LhdpF[[count]] <- Lhdp[[i]]
      lookup[[count]] <- i
      count <- count + 1
    }
  }
  ret <- list()
  ret$Lhdp <- LhdpF
  ret$data <- dataF
  ret$lookup <- lookup
  return(ret)
}
