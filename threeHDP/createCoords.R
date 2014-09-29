createCoords <- function(xcoords, ycoords, grids, i, temp_vec, temp.data, box) {
  # xcoords, ycoords, grids, i, temp_vec, temp.data
  # box : lx, rx, uy, dy
  # initialize output
  N <- length(xcoords)
  M <- length(temp.data$gaps)
  pcoords <- c()
  ncoords <- c()
  if(N == 1) {
    if (temp_vec == i) {
      xnext <- temp.data$x[M] + 0.001
      ynext <- box[3]
    } else {
      xnext <- temp.data$x[i]
      ynext <- temp.data$y[i]
    }
    ncoords <- getIntersection(xnext, ynext, box, xcoords, ycoords) 
    if( (i-1) == 1) {
      pcoords[1] <- xcoords + 0.001
      pcoords[2] <- box[3]
    } else {
      xprev <- temp.data$x[i-2]
      yprev <- temp.data$y[i-2]
      pcoords <- getIntersection(xprev, yprev, box, xcoords, ycoords) 
    }
    ret.xcoords <- c(pcoords[1], xcoords, ncoords[1])
    ret.ycoords <- c(pcoords[2], ycoords, ncoords[2])
  } else {
    if (i == M && sum(temp_vec == i)) {
      xnext <- temp.data$x[M] + 0.001
      ynext <- box[3]
    } else {
      xnext <- temp.data$x[i]
      ynext <- temp.data$y[i]
    }
    ncoords <- getIntersection(xnext, ynext, box, xcoords[N], ycoords[N])
    xprev <- temp.data$x[temp_vec[1]-1]
    yprev <- temp.data$y[temp_vec[1]-1]
    if( temp_vec[1] == 1) {
      pcoords[1] <- xcoords[1] + 0.001
      pcoords[2] <- box[3]
    } else {
      pcoords <- getIntersection(xprev, yprev, box, xcoords[1], ycoords[1])
    }
    ret.xcoords <- c(pcoords[1], xcoords, ncoords[1])
    ret.ycoords <- c(pcoords[2], ycoords, ncoords[2])
  }
  res <- list()
  res$x <- ret.xcoords
  res$y <- ret.ycoords
  return(res)
}