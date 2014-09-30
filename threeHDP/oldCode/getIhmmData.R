getIhmmData <- function(hmmdata) {
  # Test
  sapply(1:dim(hmmdata[[1]]$states)[2],  
  function(i) (hmmdata[[1]]$states[2, i] - 1)*26 + hmmdata[[1]]$states[1, i] )
}
#paths[[p.count]]$sample <- data.ihmm$sample
#paths[[p.count]]$grid <- data.ihmm$grid
#paths[[p.count]]$states <- data.ihmm$states
#paths[[p.count]]$pstates <- data.ihmm$pstates