getNbGridsF <- function(grid, cSize, Lmin, Bmin) {
  nc <- floor(Lmin / cSize) + 1
  nr <- floor(Bmin / cSize)
  lb <-  (floor(Bmin / cSize))* (floor(Lmin / cSize) + 1) + 1
  le <- (floor(Bmin / cSize) + 1) * (floor(Lmin / cSize) + 1)
  ug <- (grid - nc)
  urg <- (grid - (nc-1))
  ulg <- (grid - (nc+1))
  lg <- (grid -1)
  rg <- (grid + 1)
  drg <- (grid + (nc+1))
  dg <- grid + nc
  dlg <- grid + nc-1
  if(grid > (nc+1)) {
    res <- c(ulg, ug, urg, lg, rg, dlg, dg, drg)
  } else if(grid <= nc) {
    res <- c(lg, rg, dlg, dg, drg)
  }

  if( (grid %% nc) == 0)  {
    res <- c(ug, ulg, lg, dlg, dg)
  } else if ((grid-1) %% nc == 0) {
    res <- c(ug, urg, rg, drg, dg)
  }
  if (grid > lb) {
    res <- c(lg, ulg, ug, urg, rg)
  }
  #########################
   if (grid == lb) {
    res <- c(ug, urg, rg)
  }
  if (grid == le) {
    res <- c(ug, ulg, lg)
  }
   if (grid == 1) {
    res <- c(rg, drg, dg)
  }
  if (grid == nc) {
    res <- c(lg, dlg, dg)
  }
  return(res)
}

