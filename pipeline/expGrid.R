expGrid <- function(rtt, paths, size) {
  library("Matrix")
  library("gtools")
  source("getNbGridsF.R")
  source("resetCountsKc.R")
  source("updateCountsS.R")
  source("getTransProbsfcp.R")
  source("predHMMdtc.R")
  source("resultsScriptn.R")
  paths1 <- paths[[1]]
  gridlist <- paths[[2]]
  lookup1 <- paths[[3]]
  cSize <- size
  Lmin <- 100
  Bmin <- 80
  M <- (floor(Bmin / cSize) +1) * (floor(Lmin / cSize) + 1)
  tcountsi <- resetCountsKc(M, rtt$K, lookup1, cSize, Lmin, Bmin)
  utcounts <- updateCountsS(M, gridlist, paths1, lookup1, rtt$K, rtt$L, tcountsi, cSize)
  ptrans <- getTransProbsfcp(M, utcounts, rtt$lw, 10, rtt$K, cSize)
  startcount <- sapply(1:M, function(i) sum(sapply(1:length(paths1), function(j) paths1[[j]]$grids[1] == i)))
  
}
