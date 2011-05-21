fdr.correct.bh <- function(p.values) {
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
  m <- length(p.values)
  bonf.p <- p.values*m
  sapply(bonf.p, function(x) min(x, 1))
}
