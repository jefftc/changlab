# 080625  created

##################################################
# Extract parameters from the command line.
argv <- commandArgs()
if(length(argv) != 9)
  stop("Usage: normscript.R --vanilla <path> <annotfile> <filestem> <algorithm> <log?> <filter25?> <filter50?>")

DATA.DIR <- argv[3]
ANNOTFILE <- argv[4]
FILESTEM <- argv[5]
ALGORITHM <- argv[6]
LOG.SIGNAL <- as.numeric(argv[7])
FILTER.25 <- as.numeric(argv[8])
FILTER.50 <- as.numeric(argv[9])

MAS5.BATCH.SIZE <- 50

##DATA.DIR <- "../geo/data/datasets/GSE5460.CEL"
##ANNOTFILE <- "data/affymetrix/HG-U133_Plus_2_annot.csv.gz"
##FILESTEM <- "GSE5460"
##ALGORITHM <- "MAS5"
##LOG.SIGNAL <-TRUE
##FILTER.25 <- TRUE
##FILTER.50 <- TRUE

# Check the parameters.
if(!file.exists(DATA.DIR))
  stop("I cannot find path.")

if(!file.exists(ANNOTFILE))
  stop("I cannot find annotfile.")

ALGORITHM.U <- gsub("([a-z])", "\\U\\1", ALGORITHM, perl=TRUE)
ALGORITHM.L <- gsub("([A-Z])", "\\L\\1", ALGORITHM, perl=TRUE)
if(!length(intersect(ALGORITHM.U, c("RMA", "MAS5"))))
  stop("Algorithm should be RMA or MAS5.")

if(LOG.SIGNAL != 0 & LOG.SIGNAL != 1)
  stop("Invalid log.")

if(FILTER.25 != 0 & FILTER.25 != 1)
  stop("Invalid filter 25.")

if(FILTER.50 != 0 & FILTER.50 != 1)
  stop("Invalid filter 50.")

# Get a list of the files to analyze.
files <- list.files(DATA.DIR)
#files <- files[1:10]
if(!length(files)) stop("No files found")
filestems <- sub("\\.cel.*", "", files, perl=TRUE, ignore.case=TRUE)
#filestems <- unlist(strsplit(files, "\\..*", perl=TRUE))
fullpaths <- sapply(files, function(x) paste(DATA.DIR, "/", x, sep=""))

# Guess if the CEL files are compressed.
COMPRESSED <- length(grep("\\.gz$", files, perl=TRUE, ignore.case=TRUE)) > 0

# Generate the file names.
L2.STR <- ""
if(LOG.SIGNAL) {
  L2.STR <- ".l2"
}
FILE.TMP <- sprintf("%s.tmp.%d", FILESTEM, Sys.getpid())
FILE.NORM <- sprintf("%s.%s", FILESTEM, ALGORITHM.L)
FILE.L2 <- sprintf("%s.l2.%s", FILESTEM, ALGORITHM.L)
FILE.FILTER25 <- sprintf("%s%s.mv2525.%s", FILESTEM, L2.STR, ALGORITHM.L)
FILE.FILTER50 <- sprintf("%s%s.mv5050.%s", FILESTEM, L2.STR, ALGORITHM.L)



##################################################
# Define some helper functions.
my.read <- function(filename, header=TRUE, colClasses=NA) {
  if(length(grep("\\.gz$", filename, perl=TRUE, ignore.case=TRUE)))
    filename <- gzfile(filename)
  read.delim(filename, header=header, as.is=TRUE, comment.char="", quote="",
    colClasses=colClasses)
}

my.write <- function(X, filename, row.names=FALSE, col.names=FALSE) {
  data.out <- as.matrix(X)
  if(is.logical(row.names)) {
    if(row.names)
      row.names <- rownames(X)
    else
      row.names <- c()
  }
  if(is.logical(col.names)) {
    if(col.names)
      col.names <- colnames(X)
    else
      col.names <- c()
  }
  if(length(col.names))
    data.out <- rbind(col.names, data.out)
  if(length(row.names)) {
    if(length(col.names))
      row.names <- c("", row.names)
    data.out <- cbind(row.names, data.out)
  }
  write.table(data.out, filename, quote=FALSE, sep="\t",
    row.names=FALSE, col.names=FALSE)
}


##################################################
# Normalize the files.
library("affy")
print(sprintf("Normalizing %d files with %s.", length(fullpaths), ALGORITHM))

if(ALGORITHM.U == "RMA") {
  data.norm <- justRMA(filenames=fullpaths, compress=COMPRESSED)
} else if(ALGORITHM.U == "MAS5") {
  data.all <- c()
  for(i1 in seq(1, length(fullpaths), MAS5.BATCH.SIZE)) {
    i2 <- min(i1+MAS5.BATCH.SIZE-1, length(fullpaths))
    print(c(i1, i2))
    x <- ReadAffy(filenames=(fullpaths[i1:i2]))
    x.mas5 <- mas5(x, normalize=FALSE)
    data.all <- cbind(data.all, exprs(x.mas5))
  }
  data.norm <- new("ExpressionSet", exprs=data.all)
  data.norm <- affy.scalevalue.exprSet(data.norm)
} else {
  stop("Algorithm should be RMA or MAS5.")
}

# Write the file normalized file.
data.out <- exprs(data.norm)
data.out <- cbind(rownames(data.out), data.out)
colnames(data.out) <- c("Probe.Set.ID", filestems)
my.write(data.out, FILE.TMP, row.names=FALSE, col.names=TRUE)


##################################################
# Annotate.
print("Annotating data.")
data <- my.read(FILE.TMP, header=FALSE, colClasses="character")
sample.names <- data[1,2:(ncol(data))]
sample.names <- gsub("\"", "", sample.names)
probeset.id <- data[2:nrow(data),1]
probeset.id <- gsub("\"", "", probeset.id)
M <- data[2:nrow(data),2:ncol(data)]

# Load the annotatations.
con <- ANNOTFILE
if(length(grep("\\.gz$", ANNOTFILE, perl=TRUE, ignore.case=TRUE)))
  con <- gzfile(ANNOTFILE)
affyannot <- read.csv(con, comment.char="#", colClasses="character")

# Match the probe sets to the annotations.
I <- match(probeset.id, affyannot[["Probe.Set.ID"]])
if(any(is.na(I))) stop("mismatch")
affyannot <- affyannot[I,]
if(any(probeset.id != affyannot[["Probe.Set.ID"]])) stop("bug")

description <- affyannot[["Target.Description"]]
name <- "LocusLink"
if(!any(names(affyannot) == name))
  name <- "Entrez.Gene"
locuslink.id <- affyannot[[name]]
gene.symbol <- affyannot[["Gene.Symbol"]]

# Write the annotated data.
data.out <- cbind(probeset.id, description, locuslink.id, gene.symbol, M)
colnames(data.out) <- c(
  "Probe.Set.ID", "Description", "LocusLink", "Gene.Symbol", sample.names)
my.write(data.out, FILE.NORM, col.names=TRUE)



##################################################
# Log the data.
if(LOG.SIGNAL) {
  print("Log data")
  # Will this mess up the sample names?  (e.g. if sample is number)
  data <- my.read(FILE.NORM, header=TRUE)
  M <- as.matrix(data[,5:ncol(data)])
  M <- log2(M + 1E-10)
  data[,5:ncol(data)] <- M
  my.write(data, FILE.L2, col.names=TRUE)
  FILE.NORM <- FILE.L2
}


##################################################
# Filter data by 25%.
if(FILTER.25) {
  print("Filter 25.")
  data <- my.read(FILE.NORM, header=TRUE)
  M <- as.matrix(data[,5:ncol(data)])
  num.genes <- nrow(M)
  medians <- apply(M, 1, median)
  median.cutoff <- sort(medians)[floor(num.genes*0.25)]
  vars <- apply(M, 1, var)
  var.cutoff <- sort(vars)[floor(num.genes*0.25)]
  data <- data[medians >= median.cutoff & vars >= var.cutoff,]
  my.write(data, FILE.FILTER25, col.names=TRUE)
}


##################################################
# Filter data by 50%
if(FILTER.50) {
  print("Filtering 50.")
  data <- my.read(FILE.NORM, header=TRUE)
  M <- as.matrix(data[,5:ncol(data)])
  num.genes <- nrow(M)
  medians <- apply(M, 1, median)
  median.cutoff <- sort(medians)[floor(num.genes*0.50)]
  vars <- apply(M, 1, var)
  var.cutoff <- sort(vars)[floor(num.genes*0.50)]
  data <- data[medians >= median.cutoff & vars >= var.cutoff,]
  my.write(data, FILE.FILTER50, col.names=TRUE)
}



# Clean up the TMP file.
unlink(FILE.TMP)
