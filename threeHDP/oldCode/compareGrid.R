compareGrid <- function(d, data, hdpdata, betak, phi, clust) {
  di <- which(sapply(d, function(dk) any(dk == clust)) == TRUE)
  m <- matrix(0, nrow = 20, ncol = 26)
  for(i in 1:520) {
    if (any(di == i)) {
      x <- i %% 26
      y <- floor(i / 26) + 1
      if(x == 0) {
        x <- 26
        y <- floor(i / 26)
      }
      m[y, x] <- 1
    }
  }
  frot <- function(m) t(m)[,nrow(m):1]
  m <- frot(m)
  par(mfrow = c(2, 2), mar = c(1, 1, 1, 1), oma = c(0, 0, 2, 0), pty = "s")
  image(m, col=gray((0:255)/255))
  for (i in 2:4)
  simplePlot(hdpdata[[2]][[di[i]]], which(d[[di[i]]] == clust)[1], data, betak, di[i])
}