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
#
# find.outliers


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

find.outliers <- function(x, y, perc.init=NULL, z.cutoff=NULL, max.iter=NULL) {
  # Find outlier points.  Model the points with a linear regression
  # line, and calculate the z-score of each point.
  # Algorithm is modified from RANSAC.
  # 
  # Parameters:
  # perc.init  Percent of points to include initially.  (0.0-1.0).
  # z.cutoff   Points less than this z-score will be included in model.
  # max.iter   Maximum number of iterations before convergence.
  # 
  # Returns a list with members:
  # x      Original data (independent variable).
  # y      Original data (dependent variable).
  # y.hat  Prediction of y, given linear model.
  # coef   Coefficients of linear model.
  # rmsd   RMSD of each point (parallel to y).
  # z      Z-score of each point (parallel to y).

  if(is.null(perc.init)) perc.init <- 0.5
  if(is.null(z.cutoff)) z.cutoff <- 1
  if(is.null(max.iter)) max.iter <- 100

  if(length(x) != length(y)) stop("x and y not same length.")
  if(perc.init < 0 | perc.init > 1) stop("perc.init out of range")
  if(z.cutoff < 0.1) stop("z.cutoff too small")
  if(max.iter < 5) stop("too few iterations")

  # Initialize the model with the middle 50% of points.
  p1 <- (1.0-perc.init)/2.0
  p2 <- p1 + perc.init
  x.s <- sort(x)
  c1 <- x.s[p1*length(x)]
  c2 <- x.s[p2*length(x)]
  I.model <- (x >= c1) & (x < c2)  # which points are in the model

  max.iter <- 100
  for(i in 1:max.iter) {
    m <- glm(y[I.model] ~ x[I.model])
    y.hat <- as.numeric(cbind(1, x) %*% m$coef)
    rmsd <- sqrt((y.hat-y)**2)
    z <- rmsd / sqrt(var(rmsd))
    z[y.hat > y] <- -z[y.hat > y]

    old.model <- I.model
    I.model <- abs(z) < z.cutoff
    if(all(I.model == old.model)) break
  }

  list(x=x, y=y, y.hat=y.hat, coef=m$coef, rmsd=rmsd, z=z)
}

