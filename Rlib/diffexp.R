# find.de.genes.sam    Use SAM to find differentially expressed genes.


matrix2dataframe <- function(M) {
  x <- matrix(unlist(M), nrow(M), ncol(M))
  x <- data.frame(x, stringsAsFactors=FALSE)
  names(x) <- colnames(M)
  ow <- options("warn")  # turn off "NAs introduced by coercion"
  options(warn=-1)
  for(i in 1:ncol(x)) {
    y <- as.numeric(x[[i]])
    if(!any(is.na(y)))
      x[[i]] <- y
  }
  options(ow)
  x
}

find.de.genes.sam <- function(X, Y, DELTA, geneid=NA, genenames=NA, 
  FOLD.CHANGE=0) {
  # X must be logged.
  # Y should be the labels for two classes.  NA will be filtered out.
  # 
  # samr.plot(RESULTS$S, RESULTS$DELTA, min.foldchange=2)
  library("samr")

  I <- which(!is.na(Y))
  X <- X[,I]
  Y <- Y[I]

  Y.orig <- Y
  Y.all <- sort(unique(Y.orig))
  if(length(Y.all) != 2) stop("Y should contain exactly 2 classes.")
  Y <- rep(NA, length(Y.orig))
  Y[Y.orig == Y.all[1]] <- 1
  Y[Y.orig == Y.all[2]] <- 2
  if(any(is.na(Y))) stop("bad")

  if((length(geneid) == 1)  && is.na(geneid)) {
    ndigits <- floor(log(nrow(X), 10)) + 1
    geneid <- sprintf("GENE%0*d", ndigits, 1:nrow(X))
  }
  if((length(genenames) == 1) && is.na(genenames)) {
    genenames <- geneid
  }

  if(ncol(X) != length(Y)) stop("unaligned")
  if(nrow(X) != length(geneid)) stop("unaligned")
  if(nrow(X) != length(genenames)) stop("unaligned")

  D <- list(x=X, y=Y, logged2=TRUE, geneid=geneid, genenames=genenames)
  S <- samr(D, resp.type="Two class unpaired", nperms=100)
  DTAB <- samr.compute.delta.table(S, min.foldchange=FOLD.CHANGE)
  SIG <- samr.compute.siggenes.table(
    S, DELTA, D, DTAB, min.foldchange=FOLD.CHANGE)

  DATA <- NULL
  if((SIG$ngenes.up > 0) | (SIG$ngenes.lo > 0)) {
    SCORE <- matrix2dataframe(rbind(SIG$genes.up, SIG$genes.lo))
    # SAM messes up a few things.  Fix them.
    # The geneid and genenames are reversed.
    x <- SCORE[["Gene ID"]]
    SCORE[["Gene ID"]] <- SCORE[["Gene Name"]]
    SCORE[["Gene Name"]] <- x
    if(nrow(SCORE) >= 1) {
      # The Row is off-by-1.
      SCORE[["Row"]] <- SCORE[["Row"]]-1
    }

    x1 <- sprintf("Higher in %s", Y.all[1])
    x2 <- sprintf("Higher in %s", Y.all[2])
    Direction <- c(rep(x2, SIG$ngenes.up), rep(x1, SIG$ngenes.lo))
    x <- matrix(X[SCORE[["Row"]],], ncol=ncol(X))
    DATA <- matrix2dataframe(cbind(as.matrix(SCORE), Direction, x))
    names(DATA) <- c(names(SCORE), "Direction", Y.orig)

    if(!all(geneid[DATA[["Row"]]] == DATA[["Gene ID"]])) stop("unaligned")
  }

  list(S=S, DELTA.TAB=DTAB, SIG.TAB=SIG, DATA=DATA, DELTA=DELTA)
}

