getIhmmData1 <- function(hmmdata) {
  paths <- list()
  test <- hmmdata
  grids <- sapply(1:dim(test$states)[2],function(i) (test$states[2, i] - 1)*26 + test$states[1, i])
  paths <- extractIHMMData(test, grids)
  paths
}

# which(mapply(function(i) 100*(sum(aomw[[i]]$dt)/sum(aomw[[i]]$new.t)), 1:length(aomw)) < 15)