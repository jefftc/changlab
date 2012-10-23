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
