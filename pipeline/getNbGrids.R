getNbGrids <- function(grid) {
  ug <- (grid - 26)
  urg <- (grid - 25)
  ulg <- (grid - 27)
  lg <- (grid -1)
  rg <- (grid + 1)
  drg <- (grid + 27)
  dg <- grid + 26
  dlg <- grid + 25
  if(grid > 27) {
    res <- c(ulg, ug, urg, lg, rg, dlg, dg, drg)
  } else if(grid <= 26) {
    res <- c(lg, rg, dlg, dg, drg)
  }
  
  if( (grid %% 26) == 0)  {
    res <- c(ug, ulg, lg, dlg, dg)
  } else if ((grid-1) %% 26 == 0) {
    res <- c(ug, urg, rg, drg, dg)
  }
  if (grid > 495) {
    res <- c(lg, ulg, ug, urg, rg)
  }
  #########################
   if (grid == 495) {
    res <- c(ug, urg, rg)
  }
  if (grid == 520) {
    res <- c(ug, ulg, lg)
  }
   if (grid == 1) {
    res <- c(rg, drg, dg)
  }
  if (grid == 26) {
    res <- c(lg, dlg, dg)
  }
  return(res)
}
