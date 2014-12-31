sampleMeanjnmt <- function(data, d, j, lambda0, s, m0, J) {
  # list of class labels, each element of the list represents a cell
  K <- dim(data[[1]])[2]
  ymap <- lapply(d, function(d) which(d==j))
  initJ <- lapply(1:J,function(i) i)
  y <- mapply(function(b,c) if(length(b) != 0) c[b, ],
              ymap, data, SIMPLIFY = FALSE)
  temp.null <- NULL
  for(i in y) {if(!is.null(i)) {temp.null<- rbind(temp.null,as.matrix(i),deparse.level = 0)}}
  y <- temp.null
  if (!is.null(y)) {
    temp.m <- sum(sapply(d, function(d) sum(d == j)))
    m <- rep(0, dim(y)[2])
    for (i in 1:dim(y)[2]) {
      temp.si <- sum(y[, i])
      temp.mean <- (m0*s + (temp.si*lambda0[i]))/(temp.m*lambda0[i] + s)
      temp.var <- 1/(temp.m*lambda0[i] + s)
      m[i] <- rnorm(1, temp.mean, sqrt(temp.var))
    }
  } else {
    m <- rnorm(K, m0, sqrt(1/s))
  }
  return(m)
}
