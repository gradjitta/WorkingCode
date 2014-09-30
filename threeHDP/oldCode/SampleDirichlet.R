SampleDirichlet <- function(alpha) {
  # Generates samples from a given dirichlet distribution
  # sample ~ Dir(alpha) 
  #
  # Args:
  #   hyperparameters alpha,
  # Returns:
  #   probabilites that are dirichlet distributed with the given alphas
  out <- rep(0,length(alpha))
  norm_const <- 0
  rexp1 <- function(lambda, n) {
    u <- runif(n)
    x <- -log(u)/lambda
    return(x)
  }
  rgamma1 <- function(k, lambda) {
    y <-sum(rexp1(lambda, k)) 
    return(y)
  }
  for (i in 1:length(alpha)) {
    out[i] <- rgamma1(alpha[i], 1)
    norm_const <- norm_const + out[i]
  }
  out <- out / norm_const
  return(out)
}