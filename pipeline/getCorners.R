getCorners <- function(cSize, Lmin, Bmin) {
  nc <- floor(Lmin / cSize) + 1
  lb <-  (floor(Bmin / cSize))* (floor(Lmin / cSize) + 1) + 1
  le <- (floor(Bmin / cSize)+1) * (floor(Lmin / cSize) + 1)
  ret <- c(1,nc,lb,le)
  ret
}
