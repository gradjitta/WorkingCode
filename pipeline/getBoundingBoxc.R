getBoundingBoxc <- function(G, cSize) {
  nc <- floor(100 / cSize) + 1
  if ((G %% nc) > 0) {
    rx <- (G %% nc) * cSize
    lx <- ((G %% nc) - 1) * cSize
    uy <- -cSize * (floor(G / nc)) 
    dy <- -cSize * (ceiling(G / nc ))
    box <- c(lx, rx, uy, dy)
    } else {
      rx <- nc * cSize
      lx <- (nc-1) * cSize
      uy <- -cSize * ((G / nc) - 1) 
      dy <- -cSize * ((G / nc))
      box <- c(lx, rx, uy, dy)
      }
  box
  }

