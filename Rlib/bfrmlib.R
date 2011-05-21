# load.model
# clean.model
# load.datafile


load.model <- function(root, param.file="params.dat") {
  # Load the files in the model without any extra processing.
  # mExternalProb is relative to the original dataset.
  # param.file should be a local filename, relative to root.
  if(!file.exists(root)) stop(sprintf("Missing path %s", root))
  files <- list(
    "B.nz"=file.path(root, "mBz.txt"),
    "B"=file.path(root, "mA.txt"),
    "F"=file.path(root, "mF.txt"),
    "Psi"=file.path(root, "mPsi.txt"),
    "Tau"=file.path(root, "mTau.txt"),
    "mVariablesIn"=file.path(root, "mVariablesIn.txt"),
    "mPostPib"=file.path(root, "mPostPib.txt"),
    "mExternalProb"=file.path(root, "mExternalProb.txt")
    )
  model <- list()
  for(n in names(files)) {
    file <- files[[n]]
    if(!file.exists(file))
      next
    x <- read.delim(file, comment.char="", header=FALSE, quote="")
    # Stupid Matlab has empty columns.  Remove them.
    last.col <- x[,ncol(x)]
    if(all(is.na(last.col)))
      x <- x[,1:(ncol(x)-1)]
    model[[n]] <- as.matrix(x)
  }
  model$mVariablesIn <- as.numeric(model$mVariablesIn)
  model$Psi <- as.numeric(model$Psi)
  model$Tau <- as.numeric(model$Tau)

  if(!file.exists(param.file))
    param.file <- file.path(root, param.file)
  if(!file.exists(param.file))
    stop("missing parameter file")
  params.dat <- read.delim(param.file, 
      sep=" ", header=FALSE, as.is=TRUE, comment.char="#", quote="")

  # Set the number of designs.
  num.designs <- NA
  design.names <- c("NDesigns", "NDesignVariables")
  for(name in design.names) {
    if(any(params.dat[,1] == name))
      num.designs <- params.dat[params.dat[,1] == name, 3]
  }
  if(is.na(num.designs)) stop("Could not find design variables")
  model$NUM.DESIGNS <- as.numeric(num.designs)

  # Set the number of control variables.
  num.controls <- 0
  if(any(params.dat[,1] == "NControlVariables"))
    num.controls <- params.dat[params.dat[,1] == "NControlVariables", 3]
  model$NUM.CONTROLS <- as.numeric(num.controls)

  # Set the number of variables.
  x <- which(params.dat[,1] == "NVariables")
  if(length(x) != 1) stop("Could not find NVariables")
  model$N.VARIABLES <- as.numeric(params.dat[x,3])

  model$DataFile <- params.dat[params.dat[,1] == "DataFile", 3]

  model
}

clean.model <- function(model, FACTOR.CUTOFF=0.99) {
  # Process the model for simpler handling.
  # o Remove the vectors corresponding to the designs and controls.
  # o Ensure that mVariablesIn, even if run in non-evolutionary mode.
  # 
  # Variables:
  # o mVariablesIn   Always exists, even if run in non-evolutionary mode.
  #                  Indexes of genes, relative to original data set.
  # o mExternalProb  Probabilities, relative to original data set.
  # o factors   Gene x factor, based on FACTOR.CUTOFF.
  # o FACTOR.O  Sort the factors in decreasing size.
  # o GENE.O    Sort the genes based on factor membership.
  # 
  # Matrix variables, B, F, mPostPib, factors, mVariablesIn, etc. are
  # all sorted according to FACTOR.O and GENE.O.  FACTOR.O and GENE.O are
  # provided to help convert the original data set to the same order
  # as these variables.

  # Get rid of the design variables.
  NUM.DESIGNS <- model$NUM.DESIGNS + model$NUM.CONTROLS
  num.genes <- nrow(model$B)
  num.factors <- ncol(model$B)-NUM.DESIGNS
  num.samples <- ncol(model$F)
  B <- matrix(model$B[,(NUM.DESIGNS+1):ncol(model$B)],
    num.genes, num.factors)
  F <- matrix(model$F[(NUM.DESIGNS+1):(nrow(model$F)),],
    num.factors, num.samples)
  Tau <- model$Tau[(NUM.DESIGNS+1):length(model$Tau)]
  mPostPib <- matrix(model$mPostPib[,(NUM.DESIGNS+1):ncol(model$mPostPib)],
    num.genes, num.factors)
  if(!is.null(model$B.nz)) {
    B.nz <- matrix(model$B.nz[,(NUM.DESIGNS+1):ncol(model$B.nz)], 
      num.genes, num.factors)
  }

  # Make the factors variable.
  factors <- mPostPib
  factors[factors < FACTOR.CUTOFF] <- 0
  factors[factors > 0] <- 1

  # Set mVariablesIn if one doesn't exist.
  mVariablesIn <- model$mVariablesIn
  if(!length(mVariablesIn))
    mVariablesIn <- 1:nrow(factors)

  # order the factors based on decreasing membership
  sums <- apply(factors, 2, sum)
  FACTOR.O <- order(sums, decreasing=1)

  # Order the genes based on decreasing number of factors.  Earlier
  # factors should get much higher weights.
  weights <- 2^rev(1:ncol(factors))
  x <- factors
  x <- x[,FACTOR.O]
  y <- t(t(x)*weights)
  sums <- apply(y, 1, sum)
  GENE.O <- order(sums, decreasing=1)

  # Remove the index column in mExternalProb.
  mExternalProb <- NULL
  if(!is.null(model$mExternalProb)) {
    gene.prob <- matrix(0, model$N.VARIABLES, num.factors)
    gene.prob[model$mVariablesIn,] <- mPostPib

    if(!all(gene.prob[model$mExternalProb[,1],] == 0)) stop("overlap")
    gene.prob[model$mExternalProb[,1],] <- 
      model$mExternalProb[,2:ncol(model$mExternalProb)]
    if(any(apply(gene.prob, 1, sum)==0)) stop("0 probability")
    mExternalProb <- gene.prob
  }

  cmod <- list()
  if(!is.null(model$B.nz))
    cmod$B.nz <- matrix(B.nz[GENE.O, FACTOR.O], nrow(B.nz), ncol(B.nz))
  cmod$B <- matrix(B[GENE.O, FACTOR.O], nrow(B), ncol(B))
  cmod$mPostPib <- matrix(mPostPib[GENE.O, FACTOR.O], 
    nrow(mPostPib), ncol(mPostPib))
  cmod$factors <- matrix(factors[GENE.O, FACTOR.O], 
    nrow(factors), ncol(factors))
  cmod$F <- matrix(F[FACTOR.O,], nrow(F), ncol(F))
  cmod$Psi <- model$Psi[GENE.O]
  cmod$Tau <- Tau[FACTOR.O]
  cmod$mVariablesIn <- mVariablesIn[GENE.O]
  cmod$GENE.O <- GENE.O
  cmod$FACTOR.O <- FACTOR.O
  cmod$DataFile <- model$DataFile
  cmod$mExternalProb <- mExternalProb[,FACTOR.O]
  cmod
}

get.datafile.name <- function(cmod, FULLFILE=FALSE) {
  filename <- sprintf("../arraydata/data/%s", cmod$DataFile)
  exts <- c(".rma.gz", ".mas5.gz", ".sig.gz", ".txt.gz", ".gz")
  for(ext in exts) {
    datafile <- sub(".nude", ext, filename)
    if(file.exists(datafile)) break
  }
  if(FULLFILE) {
    filters <- c(
      ".mv5050_and", ".mv5050", ".mv0050", ".mv2525", ".zv2050", 
      ".zfill25", ".z25", ".v50", ".alle2f_minFP_good74", ".f_eindor")
    fullfile <- datafile
    for(f in filters)
      fullfile <- sub(f, "", fullfile)
    if(fullfile == datafile) 
      stop(sprintf("Could not make full filename: %s", datafile))
    datafile <- fullfile
  }
  if(!file.exists(datafile))
      stop(sprintf("Could not find file: %s", datafile))
  datafile
}

load.datafile <- function(cmod, FULLFILE=FALSE) {
  # Load the original (probably filtered) array file.  If FULLFILE is
  # TRUE, then will load the entire array file, before it has been
  # filtered.  To get just the data for the model, do:
  # data[cmod$mVariablesIn,]
  datafile <- get.datafile.name(cmod, FULLFILE=FULLFILE)
  con <- gzfile(datafile)
  data <- read.delim(con, header=TRUE, as.is=TRUE, comment.char="", quote="")
  #if(!FULLFILE)
  #  data <- data[cmod$mVariablesIn,]
  data
}

