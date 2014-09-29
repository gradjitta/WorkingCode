sampleVariancejn <- function(data, j, epi, d, mean0, J) {
  K <- dim(data[[1]])[2]
  ymap <- lapply(d, function(d) which(d==j))
  initJ <- lapply(1:J,function(i) i)
  y <- mapply(function(b,c) if(length(b) != 0) c[b, ],
              ymap, data, SIMPLIFY = FALSE)
  temp.m <- sum(sapply(d, function(d) sum(d == j)))
  temp.null <- NULL
  for(i in y) {temp.null<- rbind(temp.null,i,deparse.level = 0)}
  y <- temp.null
  if (!is.null(y)) {
    val <- rep(0, dim(y)[2])
    for (i in 1:dim(y)[2]) {
      alpha <- epi + (temp.m/2)
      temp <- sum((y[, i]-mean0[i])^2)
      beta <- epi + (temp/2)
      val[i] <- 1/rgamma(1, alpha, beta)
    }
  } else {
    val <- 1/rgamma(K, epi, epi)
  }
  return(val)
}
