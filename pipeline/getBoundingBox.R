getBoundingBox <- function(G) {
# getBoundingBoxc.R is used in the main code, Ignore this function
  if ((G %% 26) > 0) {
    rx <- (G %% 26) * 4.0
    lx <- ((G %% 26) - 1) * 4.0
    uy <- -4.0 * (floor(G / 26)) 
    dy <- -4.0 * (ceiling(G / 26 ))
    box <- c(lx, rx, uy, dy)
    } else {
      rx <- 26 * 4.0
      lx <- 25 * 4.0
      uy <- -4.0 * ((G / 26) - 1) 
      dy <- -4.0 * ((G / 26))
      box <- c(lx, rx, uy, dy)
      }
  box
  }

