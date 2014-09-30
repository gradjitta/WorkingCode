collectCounts <- function(tempCounts) {
   t1 <- tempCounts[[1]][[1]]
   t2 <- tempCounts[[1]][[2]]
   tg <- tempCounts[[1]][[3]]
   book1 <- tempCounts[[1]][[4]]
   for (p in 2:10) {
    for (i in 1:520) {
      J <- length(tempCounts[[p]][[1]][[i]])
      for (j in 1:J) {
         if (!is.null(tempCounts[[p]][[1]][[i]][[j]])) {
           t1[[i]][[j]] <- t1[[i]][[j]] + tempCounts[[p]][[1]][[i]][[j]]
           #book1[[i]][[j]] <- lapply(1:length(tempCounts[[p]][[4]][[i]][[j]]), function(jj) book1[[i]][[j]][[jj]] <- c(book1[[i]][[j]][[jj]], tempCounts[[p]][[4]][[i]][[j]][[jj]]))
         }
      }
    }
   }
   for (p in 2:10) {
      for (i in 1:520) {
         t2[[i]] <- t2[[i]] + tempCounts[[p]][[2]][[i]]
      }
  }
   for (p in 2:10) {
     tg <- tg + tempCounts[[p]][[3]]
   }
   return(list(t1,t2,tg,book1))
}
