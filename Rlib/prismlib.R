write.scatterplot <- function(filename, DATA) {
  # DATA is a list of NAME -> matrix of xy coordinates.

  #maxlen <- max(unlist(lapply(DATA, nrow)))
  data.tab <- c()
  for(i in 1:length(DATA)) {
    D <- DATA[[i]]
    x <- matrix("", nrow(D), length(DATA)+1)
    x[,1] <- D[,1]
    x[,i+1] <- D[,2]
    data.tab <- rbind(data.tab, x)
  }
  colnames(data.tab) <- c("X", names(DATA))  
  write.table(data.tab, filename, quote=FALSE, sep="\t", row.names=FALSE,
    col.names=TRUE)
}

write.boxplot <- function(filename, DATA) {
  maxlen <- max(unlist(lapply(DATA, length)))

  data.tab <- c()
  for(i in 1:length(DATA)) {
    x <- DATA[[i]]
    x <- c(x, rep("", maxlen-length(x)))
    data.tab <- cbind(data.tab, x)
  }
  colnames(data.tab) <- names(DATA)
  write.table(data.tab, filename, quote=FALSE, sep="\t", row.names=FALSE,
    col.names=TRUE)
}
