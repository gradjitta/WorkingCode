workerFunc <- function(n) {return (n^2)}
values <- 1:100
library(parallel)
numWorkers <- 8
cl <- makeCluster(numWorkers, type = "MPI")
res <- parLapply(cl, values, workerFunc)
stopClusters(cl)
mpi.exit()
print(unlist(res))
