extractIHMMData <- function(temp.data, grids) {
  # return feat, movement, grids
  temp.emit <- list()
  emit.coords <- list()
  temp.vec <- c()
  temp.coords <- matrix()
  M <- length(temp.data$x)
  new.t <- temp.data$gaps
  new.grids <- c()
  lost.t <- rep(0, M)
  tcount <- 0
  for ( i in 2:length(grids)) {
    if(grids[i-1] != grids[i]) {
      temp.vec <-c(temp.vec, i-1)
      xcoords <- temp.data$x[temp.vec]
      ycoords <- temp.data$y[temp.vec]
      N0 <- length(temp.vec)
      next.x <- temp.data$x[i]
      next.y <- temp.data$y[i]
      tlag <- temp.data$gaps[i]
      # box [lx, rx, uy, dy]
      box1 <- getBoundingBox(grids[i-1])
      box2 <- getBoundingBox(grids[i])
      new.coords1 <- createCoords(xcoords, ycoords, grids, i, temp.vec, temp.data, box1)
      coords.box1 <- getIntersection(xcoords[N0], ycoords[N0], box1, next.x, next.y)
      coords.box2 <- getIntersection(xcoords[N0], ycoords[N0], box2, next.x, next.y)
      temp.feat  <- sumsVec(new.coords1$x, new.coords1$y)
      del.t <- splitTime(c(xcoords[N0], ycoords[N0]), coords.box1, coords.box2, c(next.x, next.y), tlag)
      new.t[i-1] <- new.t[i-1] + del.t[1]
      lost.t[i] <- del.t[2]
      new.t[i]  <- del.t[3]
      tot.time <- sum(new.t[temp.vec])
      tcount <- tcount + 1
      temp.coords <- matrix(0, ncol = length(new.coords1$x), nrow = 2)
      temp.coords[1, ] <- new.coords1$x
      temp.coords[2, ] <- new.coords1$y
      emit.coords[[tcount]] <- temp.coords
      temp.emit[[tcount]] <- c(tot.time, temp.feat)
      new.grids <-c(new.grids, grids[i-1])
      temp.vec <-c()
    } else {
      temp.vec <-c(temp.vec, i-1)
    }
    if ( i == length(grids)) {
      tcount <- tcount + 1
      temp.vec <-c(temp.vec, i)
      new.grids <-c(new.grids, grids[i])
      xcoords <- temp.data$x[temp.vec]
      ycoords <- temp.data$y[temp.vec]
      box1 <- getBoundingBox(grids[i])
      new.coords1 <- createCoords(xcoords, ycoords, grids, i, temp.vec, temp.data, box1)
      temp.feat  <- sumsVec(new.coords1$x, new.coords1$y)
      temp.coords <- matrix(0, ncol = length(new.coords1$x), nrow = 2)
      temp.coords[1, ] <- new.coords1$x
      temp.coords[2, ] <- new.coords1$y
      tot.time <- sum(new.t[temp.vec])
      temp.emit[[tcount]] <- c(tot.time, temp.feat)
      emit.coords[[tcount]] <- temp.coords
      temp.vec <-c()
    }
  }
  ret.data <- list()
  ret.data$emit <- temp.emit
  ret.data$dt <- lost.t
  ret.data$grids <- new.grids
  ret.data$new.t <- new.t
  ret.data$coords <- emit.coords
  ret.data
}