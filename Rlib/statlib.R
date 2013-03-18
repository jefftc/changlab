# fdr.correct.bh
# bonferroni.correct
#
# calc.sensitivity
# calc.specificity
# calc.recall
# calc.precision
#
# logadd
# logsum
#
# normalize.quant


fdr.correct.bh <- function(p.values) {
  if(!length(p.values))
    return(NULL)
  m <- length(p.values)

  # Sort p.values in increasing order.
  O <- order(p.values)
  O.rev <- order((1:m)[O])
  p.values <- p.values[O]

  k <- 1:m
  x <- m/k * p.values
  x <- sapply(x, function(x) min(x, 1) )
  x <- sapply(1:m, function(i) min(x[i:m]) )
  x[O.rev]   # Return FDR is original order of p.values
}

bonferroni.correct <- function(p.values) {
  if(!length(p.values))
    return(NULL)
  m <- length(p.values)
  bonf.p <- p.values*m
  sapply(bonf.p, function(x) min(x, 1))
}

calc.sensitivity <- function(results) {
  # results is a vector of 1/0.  It should be sorted according to some
  # score (by decreasing score).
  cumsum(results) / sum(results)
}

calc.specificity <- function(results) {
  1 - cumsum(results==0)/sum(results==0)
}

calc.recall <- calc.sensitivity

calc.precision <- function(results, interpolate=FALSE) {
  x <- cumsum(results)/(1:length(results))
  if(interpolate) {
    #x <- sapply(1:length(x), function(i) max(x[i:length(x)]))
    for(i in (length(x)-1):1)
      x[i] <- max(x[i], x[i+1])
  }
  x
}

logadd <- function(logx, logy) {
  if(logy-logx > 100)
    return(logy)
  if(logx-logy > 100)
    return(logx)
  minxy <- min(logx, logy)
  minxy + log(exp(logx-minxy) + exp(logy-minxy))
}

logsum <- function(x) {
  x <- x[(max(x)-x) < 100]
  minx <- min(x)
  minx + log(sum(exp(x-minx)))
}

normalize.quant <- function(M) {
  #library("affy")
  library("preprocessCore")  # normalize.quantiles moved from affy library
  M.dbl <- matrix(as.double(as.matrix(M)), nrow(M), ncol(M))  # prevent crash
  N <- normalize.quantiles(M.dbl)
  colnames(N) <- colnames(M)
  N
}

