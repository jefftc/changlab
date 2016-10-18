# Microarray data.
# find.de.genes.fc
# find.de.genes.ttest            Need to import statlib.R first.
# find.de.genes.sam              Only shows significant genes.
# find.de.genes.ebayes
# find.de.genes.paired.ebayes
#
# Two-color arrays.  Comparing one color to another.
# find.de.genes.2cnoref.ebayes   Two color, only one type of sample.
#
# RNA-Seq
# find.de.genes.deseq2     Takes raw counts.
# find.de.genes.edgeR      Takes raw counts.

# Input parameters:
# X            Gene x sample matrix of logged expression values, unless
#              described otherwise.
# Y            List of classes for each sample.  Can be names.
#              Typically 2 classes only.  Will filter out NA.
# geneid       Array of gene names.
# genenames    Array of gene IDs.
# fold.change  Filter by fold change first.  This should not be logged.
#
# Each function returns a list with members:
# DATA  Table with columns:
#       Gene ID
#       Gene Name
#       p.value                 OPTIONAL
#       FDR                     OPTIONAL
#       Bonf                    OPTIONAL
#       Num Samples <Class>
#       Mean <Class>
#       Log_2 Fold Change
#       Direction
#       <Signal for Class 1>
#       <Signal for Class 2>
# Other specific to the method.




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

normalize.inputs <- function(X, Y, geneid, genenames, X.is.logged) {
  # Return a list with members:
  # n               Number of genes.
  # m               Number of samples.
  # g               Number of groups.
  # X               nxm matrix of expression values.
  # Y               m vector of 1...g
  # geneid          n vector of gene IDs.
  # genenames       n vector of gene names.
  # group.names     g vector of group names.
  # NS              g vector of num samples per group.
  # NS.not.missing  g list of n-vector of num samples for each gene.
  # X.i             g list of nxNS matrices (samples in each group)
  # MEAN            nxg matrix of the mean expression of each gene.         
  # VAR             nxg matrix of the variance of each gene.
  # lFC.group12     Fold change comparing group 1 to group 2.
  #                 NA if there are not exactly 2 groups.
  # highest         Text description of group with highest average expression.
  #                 NA if there is only 1 group.
  if(ncol(X) != length(Y)) stop("X and Y not aligned")

  # Filter out NA.
  I <- which(!is.na(Y))
  X <- X[,I]
  Y <- Y[I]

  Y.orig <- Y
  #group.names <- sort(unique(Y.orig))
  group.names <- unique(Y.orig)
  Y <- rep(NA, length(Y.orig))
  for(i in 1:length(group.names))
    Y[Y.orig == group.names[i]] <- i
  if(any(is.na(Y))) stop("missing values in Y")

  n <- nrow(X)
  m <- ncol(X)
  g <- length(group.names)
  #if(g < 2) stop("not enough groups")

  if((length(geneid) == 1)  && is.na(geneid)) {
    ndigits <- floor(log(n, 10)) + 1
    geneid <- sprintf("GENE%0*d", ndigits, 1:n)
  }
  if((length(genenames) == 1) && is.na(genenames)) {
    genenames <- geneid
  }

  # Pull out a submatrix for each class.
  cn <- colnames(X)
  if(is.null(cn))
    cn <- sprintf("SAMP%03d", 1:ncol(X))
  X.i <- list()
  for(i in 1:g) {
    x <- matrix(X[,Y==i], nrow=nrow(X))
    colnames(x) <- cn[Y==i]
    X.i[[i]] <- x
  }
  
  # Count the number of samples per group.
  NS <- rep(NA, g)
  for(i in 1:g)
    NS[i] <- sum(Y==i)
  NS.not.missing <- list()
  for(i in 1:g) {
    x <- X.i[[i]]
    num.not.missing <- apply(x, 1, function(row) sum(!is.na(row)))
    NS.not.missing[[i]] <- num.not.missing
  }

  # BUG: MEAN and VAR calculations will fail is all values are
  # missing.
  # Calculate the MEAN expression for each group.
  MEAN <- matrix(NA, n, g)
  for(i in 1:g) {
    x <- X.i[[i]]
    # will be NaN if empty matrix
    MEAN[,i] <- apply(x, 1, function(xi) mean(xi, na.rm=TRUE))
  }

  # Calculate the VARIANCE of expression for each group.
  VAR <- matrix(NA, n, g)
  for(i in 1:g) {
    x <- X.i[[i]]
    # will be NaN if empty matrix
    VAR[,i] <- apply(x, 1, function(xi) var(xi, na.rm=TRUE))
  }

  # Fold change between groups 1 and 2.
  lFC.group12 <- NA
  if(g == 2) {
    if(X.is.logged) {
      lFC.group12 <- abs(MEAN[,2] - MEAN[,1])
    } else {
      x1 <- log(X.i[[1]]+1E-300, 2)
      x2 <- log(X.i[[2]]+1E-300, 2)
      mean1 <- apply(x1, 1, mean)
      mean2 <- apply(x2, 1, mean)
      lFC.group12 <- abs(mean2 - mean1)
    }
  }

  # Figure out the highest group.
  x1 <- "Higher"
  if(g > 2)
    x1 <- "Highest"
  x2 <- rep(NA, n)
  max.i <- apply(MEAN, 1, which.max)
  for(i in 1:g)
    x2[max.i==i] <- group.names[i]
  if(any(is.na(x2))) stop("missing max values")
  highest <- NA
  if(g > 1)
    highest <- sprintf("%s in %s", x1, x2)
  # If the number is the same everywhere, then say it's SAME.
  mean.min <- apply(MEAN, 1, min)
  mean.max <- apply(MEAN, 1, max)
  highest[mean.min == mean.max] <- "SAME"

  if(nrow(X) != n) stop("unaligned 1")
  if(ncol(X) != m) stop("unaligned 2")
  if(length(Y) != m) stop("unaligned 3")
  if(length(geneid) != n) stop("unaligned 4")
  if(length(genenames) != n) stop("unaligned 5")
  if(length(group.names) != g) stop("unaligned 6")
  if(length(NS) != g) stop("unaligned 7")
  if(length(NS.not.missing) != g) stop("unaligned 8")
  for(i in 1:g)
    if(length(NS.not.missing[[i]]) != n) stop("unaligned 9")
  if(length(X.i) != g) stop("unaligned 10")
  if(nrow(MEAN) != n) stop("unaligned 11")
  if(ncol(MEAN) != g) stop("unaligned 12")
  if((g == 2) && (length(lFC.group12) != n)) stop("unaligned 13")
  if((g > 1) && (length(highest) != n)) stop("unaligned 14")

  list(n=n, m=m, g=g, X=X, Y=Y, geneid=geneid, genenames=genenames, 
    group.names=group.names, NS=NS, NS.not.missing=NS.not.missing, X.i=X.i, 
    MEAN=MEAN, VAR=VAR, lFC.group12=lFC.group12, highest=highest)
}

slice.by.genes <- function(IN, I) {
  # I is a list of indexes.
  X <- matrix(IN$X[I,], ncol=ncol(IN$X), dimnames=list(NULL, colnames(IN$X)))
  geneid <- IN$geneid[I]
  genenames <- IN$genenames[I]
  NS.not.missing <- list()
  for(i in 1:length(IN$NS.not.missing)) {
    x <- IN$NS.not.missing[[i]][I]
    NS.not.missing[[i]] <- x
  }
  X.i <- list()
  for(i in 1:length(IN$X.i)) {
    x <- matrix(IN$X.i[[i]][I,], ncol=ncol(IN$X.i[[i]]))
    colnames(x) <- colnames(IN$X.i[[i]])
    X.i[[i]] <- x
  }
  MEAN <- IN$MEAN[I,]
  VAR <- IN$VAR[I,]
  lFC.group12 <- NA
  if(IN$g == 2)
    lFC.group12 <- IN$lFC.group12[I]
  highest <- NA
  if(IN$g > 1)
    highest <- IN$highest[I]

  list(n=nrow(X), m=IN$m, g=IN$g, X=X, Y=IN$Y, 
    geneid=geneid, genenames=genenames, group.names=IN$group.names,
    NS=IN$NS, NS.not.missing=NS.not.missing,
    X.i=X.i, MEAN=MEAN, VAR=VAR, lFC.group12=lFC.group12,
    highest=highest)
}

filter.by.fold.change <- function(IN, fold.change) {
  if(fold.change < 1E-10)
    return(IN)
  I <- which(IN$lFC.group12 >= log(fold.change, 2))
  slice.by.genes(IN, I)
}

make.output.table <- function(IN, p.values, fdr, bonf, filter.p05) {
  # p.values, fdr, or bonf can be empty if not calculated by the
  # method.
  if(length(p.values))
    if(length(p.values) != IN$n) stop("bad")
  if(length(fdr))
    if(length(fdr) != IN$n) stop("bad")
  if(length(bonf))
    if(length(bonf) != IN$n) stop("bad")

  x1 <- cbind(IN$geneid, IN$genenames)
  if(length(p.values))
    x1 <- cbind(x1, p.values)
  if(length(fdr))
    x1 <- cbind(x1, fdr)
  if(length(bonf))
    x1 <- cbind(x1, bonf)
  # BUG: Should use NS.not.missing.
  x2 <- c()
  if(length(IN$NS.not.missing) != IN$g) stop("unaligned: NS.not.missing")
  for(i in 1:length(IN$NS.not.missing))
    x2 <- cbind(x2, IN$NS.not.missing[[i]])
  x3 <- cbind(IN$MEAN, IN$lFC.group12, IN$highest)
  x4 <- IN$X.i[[1]]
  if(IN$g >= 2) {
    for(i in 2:IN$g)
      x4 <- cbind(x4, IN$X.i[[i]])
  }
  DATA <- cbind(x1, x2, x3, x4)
  x1 <- c("Gene ID", "Gene Name")
  if(length(p.values))
    x1 <- c(x1, "p.value")
  if(length(fdr))
    x1 <- c(x1, "FDR")
  if(length(bonf))
    x1 <- c(x1, "Bonf")
  x2 <- c(sprintf("Num Samples %s", IN$group.names), 
    sprintf("Mean %s", IN$group.names), 
    "Log_2 Fold Change", "Direction")
  x3 <- c()
  for(i in 1:IN$g)
    x3 <- c(x3, colnames(IN$X.i[[i]]))
  x <- c(x1, x2, x3)
  # If DATA is empty, make sure it has the right number of columns.
  DATA <- matrix(DATA, ncol=length(x))
  colnames(DATA) <- x
  DATA <- matrix2dataframe(DATA)
  # Make sure Gene ID and Gene Names don't get converted to numbers
  DATA[["Gene ID"]] <- IN$geneid
  DATA[["Gene Name"]] <- IN$genenames

  if(filter.p05 && length(p.values)) {
    I <- which(p.values < 0.05)
    DATA <- DATA[I,]
  }

  DATA
}



find.de.genes.fc <- function(X, Y, geneid=NA, genenames=NA, FOLD.CHANGE=0) {
  IN <- normalize.inputs(X, Y, geneid, genenames, TRUE)
  if(IN$g != 2) stop("Y should contain exactly 2 classes.")
  # Need at least 1 sample to calculate fold change.
  if(IN$NS[1] < 1) stop("not enough samples")
  if(IN$NS[2] < 1) stop("not enough samples")

  IN <- filter.by.fold.change(IN, FOLD.CHANGE)
  DATA <- make.output.table(IN, c(), c(), c(), FALSE)
  return(list(DATA=DATA))

  # Sort by decreasing fold change.
  #O <- order(as.numeric(DATA[,2+length(Y.all)*2+1]), decreasing=TRUE)
  #DATA <- DATA[O,]
}


# requires statlib.R
find.de.genes.ttest <- function(X, Y, geneid=NA, genenames=NA, 
  FOLD.CHANGE=0, filter.p05=FALSE, NPROCS=1) {
  library(parallel)

  IN <- normalize.inputs(X, Y, geneid, genenames, TRUE)
  if(IN$g != 2) stop("Y should contain exactly 2 classes.")
  if(IN$NS[1] < 2) stop("not enough samples")
  if(IN$NS[2] < 2) stop("not enough samples")

  IN <- filter.by.fold.change(IN, FOLD.CHANGE)
  if(IN$n == 0) {
    # No genes.  Return an empty matrix.
    DATA <- make.output.table(IN, c(), c(), c(), FALSE)
    return(list(DATA=DATA))
  }

  # Calculate the p-values using a t-test.
  p.values <- rep(NA, IN$n)
  # If the variance of X.1 and X.2 are both 0, then t.test will
  # generate an error:
  # data are essentially constant.
  # If there's no variance in any of the groups, then p-value is 0.
  sumvar <- apply(IN$VAR, 1, sum)
  p.values[sumvar == 0] <- 0
  # If there's no fold change, then p-value is 1.
  p.values[IN$lFC.group12 < 1E-10] <- 1

  I <- which(is.na(p.values))
  X.1 <- IN$X.i[[1]]
  X.2 <- IN$X.i[[2]]
  # Need to filter out missing values.
  x <- mclapply(I, mc.cores=NPROCS, FUN=function(i) {
    x1 <- X.1[i,]
    x2 <- X.2[i,]
    x1 <- x1[!is.na(x1)]
    x2 <- x2[!is.na(x2)]
    # t-test requires at least 2 values in x1 and 2 in x2.
    if((length(x1) < 2) | (length(x2) < 2)) return(NA)
    t.test(x1, x2)
    })
  #x <- lapply(I, function(i) t.test(X.1[i,], X.2[i,]))
  for(i in 1:length(x)) {
    if(is.na(x[[i]])) next
    p.values[I[i]] <- x[[i]]$p.value
  }
  #p.values[I] <- unlist(lapply(x, function(x) x$p.value))

  # Ignore NA.  Might be due to missing values.
  fdr <- rep(NA, length(p.values))
  bonf <- rep(NA, length(p.values))
  #if(any(is.na(p.values))) stop("bad p.value")
  I <- !is.na(p.values)
  fdr[I] <- fdr.correct.bh(p.values[I])
  bonf[I] <- bonferroni.correct(p.values[I])

  DATA <- make.output.table(IN, p.values, fdr, bonf, filter.p05)
  list(DATA=DATA)
}

find.de.genes.sam <- function(X, Y, DELTA, geneid=NA, genenames=NA, 
  FOLD.CHANGE=0) {
  # DELTA
  # samr.plot(RESULTS$S, RESULTS$DELTA, min.foldchange=2)
  library("samr")

  IN <- normalize.inputs(X, Y, geneid, genenames, TRUE)
  if(IN$g != 2) stop("Y should contain exactly 2 classes.")
  if(IN$NS[1] < 2) stop("not enough samples")
  if(IN$NS[2] < 2) stop("not enough samples")

  D <- list(x=IN$X, y=IN$Y, logged2=TRUE, 
    geneid=IN$geneid, genenames=IN$genenames)
  S <- samr(D, resp.type="Two class unpaired", nperms=100)
  DTAB <- samr.compute.delta.table(S, min.foldchange=FOLD.CHANGE)
  SIG <- samr.compute.siggenes.table(
    S, DELTA, D, DTAB, min.foldchange=FOLD.CHANGE)

  # If there are no significant genes, return an empty matrix.
  if((SIG$ngenes.up == 0) & (SIG$ngenes.lo == 0)) {
    # No genes.  Return an empty matrix.
    IN <- slice.by.genes(IN, c())
    DATA <- make.output.table(IN, c(), c(), c(), FALSE)
    return(list(DATA=DATA, S=S, DELTA.TAB=DTAB, SIG.TAB=SIG, DELTA=DELTA))
  }

  SCORE <- matrix2dataframe(rbind(SIG$genes.up, SIG$genes.lo))
  # SAM messes up a few things.  Fix them.
  # The geneid and genenames are reversed.
  x <- SCORE[["Gene ID"]]
  SCORE[["Gene ID"]] <- SCORE[["Gene Name"]]
  SCORE[["Gene Name"]] <- x
  # The Row is off-by-1.
  if(nrow(SCORE) >= 1) {
    SCORE[["Row"]] <- SCORE[["Row"]]-1
  }

  I <- SCORE[["Row"]]
  IN <- slice.by.genes(IN, I)

  #x1 <- sprintf("Higher in %s", Y.all[1])
  #x2 <- sprintf("Higher in %s", Y.all[2])
  #Direction <- c(rep(x2, SIG$ngenes.up), rep(x1, SIG$ngenes.lo))
  #x <- matrix(X[SCORE[["Row"]],], ncol=ncol(X))
  #DATA <- matrix2dataframe(cbind(as.matrix(SCORE), Direction, x))
  #names(DATA) <- c(names(SCORE), "Direction", Y.orig)

  # Get the FDR values.
  I <- which(names(SCORE) == "q-value(%)")
  if(length(I) != 1) stop("missing q-value(%)")
  fdr <- SCORE[,I]/100
  #nlfdr <- -log(DATA[,I]/100+1E-25, 10)
  #DATA <- cbind(DATA[,1:3], nlfdr, DATA[,4:ncol(DATA)])
  #names(DATA)[4] <- "NL10 FDR"

  # Convert "Fold Change" to log_2.
  #I <- which(names(DATA) == "Fold Change")
  #if(length(I) != 1) stop("missing Fold Change")
  #lfc <- abs(log(DATA[,I], 2))
  #DATA[,I] <- lfc
  #names(DATA)[I] <- "Log_2 Fold Change"

  #if(!all(geneid[DATA[["Row"]]] == DATA[["Gene ID"]])) stop("unaligned")

  #list(DATA=DATA, S=S, DELTA.TAB=DTAB, SIG.TAB=SIG, DELTA=DELTA)

  DATA <- make.output.table(IN, c(), fdr, c(), FALSE)
  list(DATA=DATA, S=S, DELTA.TAB=DTAB, SIG.TAB=SIG, DELTA=DELTA, SCORE=SCORE)
}


find.de.genes.ebayes <- function(X, Y, geneid=NA, genenames=NA, 
  FOLD.CHANGE=0, filter.p05=FALSE) {
  # Can generate Warning message:
  # Zero sample variances detected, have been offset 
  library(limma)

  IN <- normalize.inputs(X, Y, geneid, genenames, TRUE)
  if(IN$g != 2) stop("Y should contain exactly 2 classes.")
  if(IN$NS[1] < 2) stop("not enough samples")
  if(IN$NS[2] < 2) stop("not enough samples")

  # Use limma to calculate statistics.
  GROUP1 <- as.numeric(IN$Y == 1)
  GROUP2 <- as.numeric(IN$Y == 2)
  design <- cbind(GROUP1=GROUP1, GROUP2=GROUP2)
  # allows for missing values!
  fit <- lmFit(IN$X, design=design)
  fit2 <- contrasts.fit(fit, c(-1, 1))
  fit2 <- eBayes(fit2)
  lfc <- 0
  if(FOLD.CHANGE > 0)
    lfc <- log(FOLD.CHANGE, 2)
  # p.value cutoff is for adjusted p-values.  We want non-adjusted p-values.
  #TOP <- topTable(fit2, adjust="fdr", lfc=log(FOLD.CHANGE, 2), p.value=0.05,
  #  number=nrow(fit2))
  TOP <- topTable(fit2, adjust="fdr", number=nrow(fit2), lfc=lfc)
  #write.table(TOP, "out.dat", col.names=TRUE, sep="\t", quote=FALSE)
  #if(!all.genes) {
  #  TOP <- TOP[TOP[["P.Value"]] < 0.05,]
  #}

  if(!nrow(TOP)) {
    # No genes.  Return an empty matrix.
    IN <- slice.by.genes(IN, c())
    DATA <- make.output.table(IN, c(), c(), c(), FALSE)
    return(list(DATA=DATA))
  }

  p.value <- TOP[["P.Value"]]
  fdr <- TOP[["adj.P.Val"]]
  #bonf <- bonferroni.correct(p.value)

  I <- as.numeric(rownames(TOP))
  IN <- slice.by.genes(IN, I)

  DATA <- make.output.table(IN, p.value, fdr, c(), filter.p05)
  list(DATA=DATA, TOP=TOP)
}


# Assume that the pairing is in order.
find.de.genes.paired.ebayes <- function(X, Y, geneid=NA, genenames=NA, 
  FOLD.CHANGE=0, filter.p05=FALSE) {
  library(limma)

  IN <- normalize.inputs(X, Y, geneid, genenames, TRUE)
  if(IN$g != 2) stop("Y should contain exactly 2 classes.")
  if(IN$NS[1] < 2) stop("not enough samples")
  if(IN$NS[2] < 2) stop("not enough samples")

  # Use limma to calculate statistics.
  GROUP1 <- as.numeric(IN$Y == 1)
  GROUP2 <- as.numeric(IN$Y == 2)

  # For paired analysis, number of samples in the classes should be
  # the same.
  if(sum(GROUP1 == 1) != sum(GROUP2 == 1)) stop("unequal number of samples")

  # Assume samples are in paired order.
  num.pairs <- sum(GROUP1)
  PAIRS <- rep(0, length(GROUP1))
  PAIRS[GROUP1 == 1] <- 1:num.pairs
  PAIRS[GROUP2 == 1] <- 1:num.pairs
  if(any(PAIRS == 0)) stop("bad")
  PAIRS <- factor(PAIRS)
  Y.fact <- factor(IN$group.names[IN$Y], levels=IN$group.names)
  design <- model.matrix(~PAIRS+Y.fact)
  # Get message:
  # Removing intercept from test coefficients
  fit <- lmFit(X, design=design)
  fit2 <- eBayes(fit)
  lfc <- 0
  if(FOLD.CHANGE > 0)
    lfc <- log(FOLD.CHANGE, 2)
  # p.value cutoff is for adjusted p-values.  We want non-adjusted p-values.
  #TOP <- topTable(fit2, adjust="fdr", lfc=log(FOLD.CHANGE, 2), p.value=0.05,
  #  number=nrow(fit2))
  TOP <- topTable(fit2, adjust="fdr", number=nrow(fit2), lfc=lfc)
  #write.table(TOP, "out.dat", col.names=TRUE, sep="\t", quote=FALSE)
  #if(!all.genes) {
  #  TOP <- TOP[TOP[["P.Value"]] < 0.05,]
  #}
  if(!nrow(TOP)) {
    # No genes.  Return an empty matrix.
    IN <- slice.by.genes(IN, c())
    DATA <- make.output.table(IN, c(), c(), c(), FALSE)
    return(list(DATA=DATA))
  }

  I <- as.numeric(rownames(TOP))
  IN <- slice.by.genes(IN, I)

  p.value <- TOP[["P.Value"]]
  fdr <- TOP[["adj.P.Val"]]
  #bonf <- bonferroni.correct(p.value)

  DATA <- make.output.table(IN, p.value, fdr, c(), filter.p05)
  list(DATA=DATA, fit=fit2, TOP=TOP)
}


# 2-color no reference
find.de.genes.2cnoref.ebayes <- function(X, geneid=NA, genenames=NA, 
  FOLD.CHANGE=0, filter.p05=FALSE) {
  library(limma)

  if(length(FOLD.CHANGE) != 1) stop("bad arguments")

  Y <- rep(1, ncol(X))  # All samples are the same class.
  IN <- normalize.inputs(X, Y, geneid, genenames, TRUE)
  if(IN$g != 1) stop("Should be only 1 class.")

  fit <- lmFit(IN$X)
  fit2 <- eBayes(fit)
  lfc <- 0
  if(FOLD.CHANGE > 0)
    lfc <- log(FOLD.CHANGE, 2)
  TOP <- topTable(fit2, adjust="fdr", number=nrow(fit2), lfc=lfc)

  if(!nrow(TOP)) {
    # No genes.  Return an empty matrix.
    IN <- slice.by.genes(IN, c())
    DATA <- make.output.table(IN, c(), c(), c(), FALSE)
    return(list(DATA=DATA, fit=fit2, TOP=TOP))
  }

  I <- as.numeric(rownames(TOP))
  IN <- slice.by.genes(IN, I)

  p.value <- TOP[["P.Value"]]
  fdr <- TOP[["adj.P.Val"]]
  #bonf <- bonferroni.correct(p.value)

  DATA <- make.output.table(IN, p.value, fdr, c(), filter.p05)
  list(DATA=DATA, fit=fit2, TOP=TOP)
}


find.de.genes.deseq2 <- function(X, Y, geneid=NA, genenames=NA, 
  FOLD.CHANGE=0, filter.p05=FALSE, NPROCS=1) {
  # X should be unnormalized counts per gene
  # summary(res)
  # plotMA(res, main="DESeq2", ylim=c(-2,2))

  library(DESeq2)
  if(NPROCS > 1) {
    library(BiocParallel)
    register(MulticoreParam(NPROCS))
  }

  IN <- normalize.inputs(X, Y, geneid, genenames, FALSE)
  if(IN$g != 2) stop("Y should contain exactly 2 classes.")
  if(IN$NS[1] < 2) stop("not enough samples")
  if(IN$NS[2] < 2) stop("not enough samples")

  IN <- filter.by.fold.change(IN, FOLD.CHANGE)
  if(IN$n == 0) {
    # No genes.  Return an empty matrix.
    DATA <- make.output.table(IN, c(), c(), c(), FALSE)
    return(list(DATA=DATA))
  }

  count.data <- IN$X
  condition <- IN$group.names[IN$Y]
  col.data <- data.frame(condition=condition)
  dds <- DESeqDataSetFromMatrix(
    countData=count.data, colData=col.data, design=~condition)  
  dds.2 <- DESeq(dds)
  # DESeq does shrinkage, which throws off fold change.  Also look at
  # the fold change values without shrinkage.
  # Contains the MLE estimate of fold change.
  # Include a column: lfcMLE
  res <- results(dds.2, addMLE=TRUE)
  # Describes columns in res
  # mcols(res)$description
  #print(mcols(res)$description)
  #write.table(res, "out2.dat", col.names=TRUE)

  p.value <- res[["pvalue"]]
  fdr <- res[["padj"]]
  #bonf <- bonferroni.correct(p.value)
  fold.change <- res[["log2FoldChange"]]
  fold.change.mle <- res[["lfcMLE"]]

  DATA <- make.output.table(IN, p.value, fdr, c(), filter.p05)
  DATA[["Log_2 Fold Change"]] <- abs(fold.change)
  # Insert the MLE fold change.
  i <- which(names(DATA) == "Log_2 Fold Change")
  if(length(i) == 0)
    stop("could not find fold change")
  DATA <- cbind(DATA[,1:i], abs(fold.change.mle), DATA[,(i+1):ncol(DATA)])
  names(DATA)[i+1] <- "Log_2 Fold Change (MLE)"
  list(DATA=DATA, dds=dds.2, res=res)
  #list(DATA=DATA, dds=dds.2, res=res, dds.noshrink=dds.2.noshrink, 
  #  res.noshrink=res.noshrink)
}


find.de.genes.edgeR <- function(X, Y, geneid=NA, genenames=NA, 
  FOLD.CHANGE=0, filter.p05=FALSE, tagwise.dispersion=TRUE) {
  # X should be unnormalized counts per gene
  # Use tagwise dispersion if there are enough replicates to estimate
  # the dispersion for individual genes (or try them both).

  library(edgeR)

  IN <- normalize.inputs(X, Y, geneid, genenames, FALSE)
  if(IN$g != 2) stop("Y should contain exactly 2 classes.")
  if(IN$NS[1] < 2) stop("not enough samples")
  if(IN$NS[2] < 2) stop("not enough samples")

  IN <- filter.by.fold.change(IN, FOLD.CHANGE)
  if(IN$n == 0) {
    # No genes.  Return an empty matrix.
    DATA <- make.output.table(IN, c(), c(), c(), FALSE)
    return(list(DATA=DATA))
  }

  group <- factor(IN$Y)
  dge <- DGEList(counts=IN$X, group=group)
  # Normalize by TMM.
  dge <- calcNormFactors(dge)
  dge <- estimateCommonDisp(dge)
  if(tagwise.dispersion)
    # OPTIONAL: Do tagwise, rather than common dispersion.
    dge <- estimateTagwiseDisp(dge)
  et <- exactTest(dge)
  top <- topTags(et, n=Inf)

  O <- order(as.numeric(rownames(top$table)))
  top.table <- top$table[O,]

  if(IN$n != nrow(top.table)) stop("unaligned")
  p.value <- top.table$PValue
  fdr <- top.table$FDR
  #bonf <- bonferroni.correct(p.value)

  DATA <- make.output.table(IN, p.value, fdr, c(), filter.p05)
  DATA[["Log_2 Fold Change"]] <- abs(top.table$logFC)
  list(DATA=DATA, dge=dge, top=top)
}
