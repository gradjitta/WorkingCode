getDist <- function(g1,g2,gn, cSize) {
   source("getBoundingBoxc.R")
   ln <- getBoundingBoxc(gn, cSize)
   l1 <- getBoundingBoxc(g1, cSize)
   l2 <- getBoundingBoxc(g2, cSize)
   dist1 <- sqrt( ((ln-l1)[1])^2 +  ((ln-l1)[3] )^2)
   dist2 <- sqrt( ((ln-l2)[1])^2 +  ((ln-l2)[3] )^2)
   ret <- c(dist1, dist2)
   ret
}
