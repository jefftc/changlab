library("SparseM")

sparse2list <- function(m, pad=0) {
  x <- list(
    ra=m@ra, ja=m@ja, ia=m@ia,
    dimension=attr(m, "dimension"),
    klass=class(m))

  if(pad) {
    # Need to pad out each subcomponent of the list so that it can be
    # saved to a file with write.table.
    l <- max(sapply(x, length))
    for(i in 1:length(x))
      x[[i]] <- c(x[[i]], rep(NA, l-length(x[[i]])))
  }
  x
}

list2sparse <- function(l) {
  # Unpad the elements of the list.
  x <- list()
  for(i in 1:length(l))
    x[[i]] <- l[[i]][!is.na(l[[i]])]
  names(x) <- names(l)
  l <- x

  # Make sure "ia" is monotonically increasing.
  if(!all(l$ia == sort(l$ia))) {
    stop("ERROR: ia is not monotonically increasing.")
  }

  x <- data.frame()
  class(x) <- as.character(l$klass)
  x@ra <- l$ra; x@ja <- l$ja; x@ia <- l$ia
  x@dimension <- l$dimension
  # Takes a long time...
  x <- x[,]   # Make it into a real sparse matrix (otherwise problems occur).
  x
}


