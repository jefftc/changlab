# write.scatterplot
# write.boxplot

write.scatterplot <- function(filename, DATA, rownames=NA) {
  # DATA is a list of NAME -> matrix of xy coordinates.
  # rownames is a list of NAME -> list of names (same length as nrow(matrix)).
  # One x, lots of Y.

  #maxlen <- max(unlist(lapply(DATA, nrow)))
  data.tab <- c()
  rn <- c()
  R <- NA
  for(i in 1:length(DATA)) {
    D <- DATA[[i]]
    if(!all(is.na(rownames))) {
      R <- rownames[[names(DATA)[i]]]
      if(is.null(R)) stop("missing name")
      if(length(R) != nrow(D)) stop("different lengths")
    }
    x <- matrix("", nrow(D), length(DATA)+1)
    x[,1] <- D[,1]
    x[,i+1] <- D[,2]
    data.tab <- rbind(data.tab, x)
    rn <- c(rn, R)
  }
  colnames(data.tab) <- c("X", names(DATA))
  if(!all(is.na(rownames))) {
    data.tab <- cbind(rn, data.tab)
    colnames(data.tab) <- c("Title", colnames(data.tab)[2:ncol(data.tab)])
  }
  write.table(data.tab, filename, quote=FALSE, sep="\t", row.names=FALSE,
    col.names=TRUE)
}

write.boxplot <- function(filename, DATA) {
  # DATA is list of name -> vector of values.
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
