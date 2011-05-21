# MOVED TO plotfns.R.
#normalize.one.mv <- function(x, M=0, V=1) {
#  # Normalize a list of numbers so the mean is M and variance is V.
#  V.0 <- var(x)
#  M.0 <- mean(x)
#  (x-M.0)*sqrt(V/V.0) + M
#}
#
#normalize.mv <- function(X, M=0, V=1) {
#  t(apply(X, 1, function(x) normalize.one.mv(x, M, V)))
#}

normalize.quant <- function(M) {
  #library("affy")
  library("preprocessCore")  # normalize.quantiles moved from affy library
  M.dbl <- matrix(as.double(as.matrix(M)), nrow(M), ncol(M))  # prevent crash
  N <- normalize.quantiles(M.dbl)
  colnames(N) <- colnames(M)
  N
}

load.and.merge.datasets <- function(file1, file2, map.file) {
  data.1 <- read.delim(file1,
    header=TRUE, as.is=TRUE, comment.char="", quote="")
  data.2 <- read.delim(file2,
      header=TRUE, as.is=TRUE, comment.char="", quote="")
  data.map <- read.delim(map.file, 
    header=FALSE, as.is=TRUE, comment.char="", quote="")
  data.1 <- data.1[data.map[,1],]
  data.2 <- data.2[data.map[,2],]
  list(data.1, data.2)
}
