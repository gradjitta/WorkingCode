getClusters <- function(d, data, hdpdata, betak, phi, clust) {
  mapC <- which(sapply(d, function(dk) any(dk == clust)) == TRUE)
  N0 <- length(mapC)
  print(betak[clust])
  print(phi[clust, ])
  grid.cells <- sapply(hdpdata[[2]], function(d) d)
  print(grid.cells)
  par(mfrow = c(5, 5), mar = c(1, 1, 1, 1), oma = c(0, 0, 2, 0), pty = "s")
  for (i in mapC) {
    simplePlot(hdpdata[[2]][[i]], which(d[[i]] == clust)[1], data, betak, i)
  }
}