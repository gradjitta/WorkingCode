getDist <- function(g1,g2,gn) {
   source("getBoundingBox.R")
   ln <- getBoundingBox(gn)
   l1 <- getBoundingBox(g1)
   l2 <- getBoundingBox(g2)
   dist1 <- sqrt( ((ln-l1)[1])^2 +  ((ln-l1)[3] )^2)
   dist2 <- sqrt( ((ln-l2)[1])^2 +  ((ln-l2)[3] )^2)
   ret <- c(dist1, dist2)
   ret
}
