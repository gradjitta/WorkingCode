getUsedClusters <- function(M, tcount, K) {
  final <- rep(0, K)
  for(i in 1:M) {
    cNULL <- sum(sapply(1:length(tcount[[i]]), function(pl) is.null(tcount[[i]][[pl]])))
    if (cNULL  == 0) {
      for (j in 1: length(tcount[[i]])) {
        ztimes <- K-dim(tcount[[i]][[j]])[2]
        if (length(ztimes) != 0) {
          padding <- rep(0, K-dim(tcount[[i]][[j]])[2])
        } else {
          padding <- c()
        }
        final <- final + c(colSums(tcount[[i]][[j]]), padding)
      }
    }
  }
  final
}
