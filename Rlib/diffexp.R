# find.de.genes.sam    Use SAM to find differentially expressed genes.
# find.de.genes.ttest
# find.de.genes.ebayes
# find.de.genes.2cnoref.ebayes


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
  # FOLD.CHANGE should not be logged.
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

    # Convert "Fold Change" to log_2.
    I <- which(names(DATA) == "Fold Change")
    if(length(I) != 1) stop("missing Fold Change")
    lfc <- abs(log(DATA[,I], 2))
    DATA[,I] <- lfc
    names(DATA)[I] <- "Log_2 Fold Change"

    if(!all(geneid[DATA[["Row"]]] == DATA[["Gene ID"]])) stop("unaligned")
  }

  list(S=S, DELTA.TAB=DTAB, SIG.TAB=SIG, DATA=DATA, DELTA=DELTA)
}


# requires statlib.R
find.de.genes.ttest <- function(X, Y, geneid=NA, genenames=NA, 
  FOLD.CHANGE=0, NPROCS=1) {
  library(multicore)
  # X must be logged.
  # Y should be the labels for two classes.  NA will be filtered out.
  # FOLD.CHANGE should not be logged.

  if(ncol(X) != length(Y)) stop("X and Y not aligned")

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
  if(sum(Y == 1, na.rm=TRUE) <= 1) stop("not enough samples")
  if(sum(Y == 2, na.rm=TRUE) <= 1) stop("not enough samples")

  if((length(geneid) == 1) && is.na(geneid)) {
    ndigits <- floor(log(nrow(X), 10)) + 1
    geneid <- sprintf("GENE%0*d", ndigits, 1:nrow(X))
  }
  if((length(genenames) == 1) && is.na(genenames)) {
    genenames <- geneid
  }
  
  if(ncol(X) != length(Y)) stop("unaligned")
  if(nrow(X) != length(geneid)) stop("unaligned")
  if(nrow(X) != length(genenames)) stop("unaligned")

  X.1 <- X[,Y == 1]
  X.2 <- X[,Y == 2]
  if(is.null(colnames(X.1)))
    colnames(X.1) <- sprintf("S1_%03d", 1:ncol(X.1))
  if(is.null(colnames(X.2)))
    colnames(X.2) <- sprintf("S2_%03d", 1:ncol(X.2))


  # Find the genes that match the fold change criteria.
  med.1 <- apply(X.1, 1, mean)
  med.2 <- apply(X.2, 1, mean)
  diff <- abs(med.2 - med.1)
  I <- diff >= log(FOLD.CHANGE, 2)
  X <- matrix(X[I,], ncol=ncol(X), dimnames=list(NULL, colnames(X)))
  X.1 <- matrix(X.1[I,], ncol=ncol(X.1), dimnames=list(NULL, colnames(X.1)))
  X.2 <- matrix(X.2[I,], ncol=ncol(X.2), dimnames=list(NULL, colnames(X.2)))
  geneid <- geneid[I]
  genenames <- genenames[I]
  med.1 <- med.1[I]
  med.2 <- med.2[I]
  diff <- diff[I]

  # BUG: will break if p-value is NA.
  p.values <- c()
  nl10p <- c()
  if(nrow(X.1) >= 1) {
    #x <- lapply(1:nrow(X.1), function(i) t.test(X.1[i,], X.2[i,]))
    x <- mclapply(1:nrow(X.1), mc.cores=NPROCS, FUN=function(i) 
      t.test(X.1[i,], X.2[i,]))
    p.values <- unlist(lapply(x, function(x) x$p.value))
  }
  fdr <- fdr.correct.bh(p.values)
  bonf <- bonferroni.correct(p.values)

  I <- which(p.values < 0.05)
  X <- matrix(X[I,], ncol=ncol(X), dimnames=list(NULL, colnames(X)))
  X.1 <- matrix(X.1[I,], ncol=ncol(X.1), dimnames=list(NULL, colnames(X.1)))
  X.2 <- matrix(X.2[I,], ncol=ncol(X.2), dimnames=list(NULL, colnames(X.2)))
  geneid <- geneid[I]
  genenames <- genenames[I]
  med.1 <- med.1[I]
  med.2 <- med.2[I]
  diff <- diff[I]
  p.values <- p.values[I]
  fdr <- fdr[I]
  bonf <- bonf[I]

  nl10p <- c()
  nl10fdr <- c()
  nl10bonf <- c()
  if(length(p.values)) {
    nl10p <- -log(p.values, 10)
    nl10fdr <- -log(fdr, 10)  
    nl10bonf <- -log(bonf, 10)
  }

  direction <- rep(sprintf("Higher in %s", Y.all[1]), length(med.1))
  direction[med.2 > med.1] <- sprintf("Higher in %s", Y.all[2])
  DATA <- cbind(geneid, genenames, nl10p, nl10fdr, nl10bonf, diff, 
    direction, X.1, X.2)
  colnames(DATA) <- c("Gene ID", "Gene Name", "NL10P", 
    "NL10 FDR", "NL10 Bonf", "Log_2 Fold Change", "Direction", 
    colnames(X.1), colnames(X.2))
  DATA <- matrix2dataframe(DATA)
  list(DATA=DATA)
}

find.de.genes.ebayes <- function(X, Y, geneid=NA, genenames=NA, 
  FOLD.CHANGE=0) {
  library(limma)
  # X must be logged.
  # Y should be the labels for two classes.  NA will be filtered out.

  if(ncol(X) != length(Y)) stop("X and Y not aligned")

  I <- which(!is.na(Y))
  X <- X[,I]
  Y <- Y[I]

  Y.all <- sort(unique(Y))
  if(length(Y.all) != 2) stop("Y should contain exactly 2 classes.")

  GROUP1 <- rep(0, length(Y))
  GROUP2 <- rep(0, length(Y))
  GROUP1[Y == Y.all[1]] <- 1
  GROUP2[Y == Y.all[2]] <- 1
  if(sum(GROUP1 == 1, na.rm=TRUE) <= 1) stop("not enough samples")
  if(sum(GROUP2 == 1, na.rm=TRUE) <= 1) stop("not enough samples")

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

  design <- cbind(GROUP1=GROUP1, GROUP2=GROUP2)
  fit <- lmFit(X, design=design)
  fit2 <- contrasts.fit(fit, c(-1, 1))
  fit2 <- eBayes(fit2)
  TOP <- topTable(fit2, adjust="fdr", lfc=log(FOLD.CHANGE, 2), p.value=0.05,
    number=nrow(fit2))
  if(!nrow(TOP)) {
    DATA <- matrix(NA, 0, 7+ncol(X))
    colnames(DATA) <- c("Gene ID", "Gene Name", "NL10P", 
      "NL10 FDR", "NL10 Bonf", "Log_2 Fold Change", "Direction", 
      colnames(X))
    DATA <- matrix2dataframe(DATA)
    return(list(DATA=DATA, fit=fit2, TOP=TOP))
  }

  nl10p <- -log(TOP[["P.Value"]], 10)
  nl10fdr <- -log(TOP[["adj.P.Val"]], 10)
  bonf <- bonferroni.correct(TOP[["P.Value"]])
  nl10bonf <- -log(bonf, 10)
  diff <- abs(TOP[["logFC"]])
  I <- as.numeric(rownames(TOP))
  X.1 <- X[I, Y == Y.all[1]]
  X.2 <- X[I, Y == Y.all[2]]
  if(is.null(colnames(X.1)))
    colnames(X.1) <- sprintf("S1_%03d", 1:ncol(X.1))
  if(is.null(colnames(X.2)))
    colnames(X.2) <- sprintf("S2_%03d", 1:ncol(X.2))

  direction <- rep(sprintf("Higher in %s", Y.all[1]), nrow(TOP))
  direction[TOP[["logFC"]] > 0] <- sprintf("Higher in %s", Y.all[2])
  DATA <- cbind(geneid[I], genenames[I], nl10p, nl10fdr, nl10bonf, diff, 
    direction, X.1, X.2)
  colnames(DATA) <- c("Gene ID", "Gene Name", "NL10P", 
    "NL10 FDR", "NL10 Bonf", "Log_2 Fold Change", "Direction", 
    colnames(X.1), colnames(X.2))
  DATA <- matrix2dataframe(DATA)
  list(DATA=DATA, fit=fit2, TOP=TOP)
}



find.de.genes.2cnoref.ebayes <- function(X, geneid=NA, genenames=NA, 
  FOLD.CHANGE=0) {
  library(limma)
  # X must be logged.

  if((length(geneid) == 1)  && is.na(geneid)) {
    ndigits <- floor(log(nrow(X), 10)) + 1
    geneid <- sprintf("GENE%0*d", ndigits, 1:nrow(X))
  }
  if((length(genenames) == 1) && is.na(genenames)) {
    genenames <- geneid
  }
  if(nrow(X) != length(geneid)) stop("unaligned")
  if(nrow(X) != length(genenames)) stop("unaligned")

  fit <- lmFit(X)
  fit <- eBayes(fit)
  TOP <- topTable(fit, adjust="fdr", lfc=log(FOLD.CHANGE, 2), p.value=0.05,
    number=nrow(fit))
  if(!nrow(TOP)) {
    DATA <- matrix(NA, 0, 7+ncol(X))
    colnames(DATA) <- c("Gene ID", "Gene Name", "NL10P", 
      "NL10 FDR", "NL10 Bonf", "Log_2 Fold Change", "Direction", 
      colnames(X))
    DATA <- matrix2dataframe(DATA)
    return(list(DATA=DATA, fit=fit, TOP=TOP))
  }

  nl10p <- -log(TOP[["P.Value"]], 10)
  nl10fdr <- -log(TOP[["adj.P.Val"]], 10)
  bonf <- bonferroni.correct(TOP[["P.Value"]])
  nl10bonf <- -log(bonf, 10)
  diff <- abs(TOP[["logFC"]])
  direction <- rep("High", nrow(TOP))
  direction[TOP[["logFC"]] < 0] <- "Low"

  I <- as.numeric(rownames(TOP))
  DATA <- cbind(geneid[I], genenames[I], nl10p, nl10fdr, nl10bonf, diff, 
    direction, X[I,])
  colnames(DATA) <- c("Gene ID", "Gene Name", "NL10P", 
    "NL10 FDR", "NL10 Bonf", "Log_2 Fold Change", "Direction", 
    colnames(X))
  DATA <- matrix2dataframe(DATA)
  list(DATA=DATA, fit=fit, TOP=TOP)
}
