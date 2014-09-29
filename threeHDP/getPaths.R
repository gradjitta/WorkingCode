getPaths <- function(data) {
  paths <- list()
  temp.count <- 0
  for (day in 1:31) {
    print("day")
    print(day)
    if (any(data[, 1] == day)) {
      dota <- genHmmData(data, day)
      p.count <- 0
      for (j in 1:length(dota)) {
        print(j)
        if ( length(dota[[j]]$ids) > 5 ) {
          p.count <- p.count + 1
          temp.paths <- getIhmmData1(dota[[j]])
          paths[[p.count + temp.count]] <- list()
          paths[[p.count + temp.count]]$emit <- temp.paths$emit
          paths[[p.count + temp.count]]$grids <- temp.paths$grids
          #paths[[p.count + temp.count]]$coords <- temp.paths[[i]]$coords
          #paths[[p.count + temp.count]]$t <- temp.paths[[i]]$new.t
          #paths[[p.count + temp.count]]$dt <- temp.paths[[i]]$dt
        }
      }
      temp.count <-temp.count + p.count
    }
  }
  return(paths)
}