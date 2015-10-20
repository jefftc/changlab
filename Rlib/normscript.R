# 080625  created

##################################################
# Extract parameters from the command line.
argv <- commandArgs()
if(length(argv) != 7)
  stop("Usage: normscript.R --vanilla <path> <annotfile> <filestem> <algorithm> <oligo>")

DATA.DIR <- argv[3]
ANNOTFILE <- argv[4]
FILESTEM <- argv[5]
ALGORITHM <- argv[6]
IS.OLIGO <- as.numeric(argv[7])

MAS5.BATCH.SIZE <- 100

##DATA.DIR <- "../geo/data/datasets/GSE5460.CEL"
##ANNOTFILE <- "data/affymetrix/HG-U133_Plus_2_annot.csv.gz"
##FILESTEM <- "GSE5460"
##ALGORITHM <- "MAS5"

# Check the parameters.
if(!file.exists(DATA.DIR))
  stop("I cannot find path.")

if(!file.exists(ANNOTFILE))
  stop("I cannot find annotfile.")

ALGORITHM.U <- gsub("([a-z])", "\\U\\1", ALGORITHM, perl=TRUE)
ALGORITHM.L <- gsub("([A-Z])", "\\L\\1", ALGORITHM, perl=TRUE)
if(!length(intersect(ALGORITHM.U, c("RMA", "MAS5"))))
  stop("Algorithm should be RMA or MAS5.")

if(IS.OLIGO != 0 & IS.OLIGO != 1)
  stop("Invalid value for oligo.  Should be 0 or 1.")

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
FILE.TMP <- sprintf("%s.tmp.%d", FILESTEM, Sys.getpid())
FILE.NORM <- sprintf("%s.%s", FILESTEM, ALGORITHM.L)



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
print(sprintf("Normalizing %d files with %s.", length(fullpaths), ALGORITHM))

if(ALGORITHM.U == "RMA" & !IS.OLIGO) {
  library(affy)
  data.norm <- justRMA(filenames=fullpaths, compress=COMPRESSED)
} else if(ALGORITHM.U == "RMA" & IS.OLIGO) {
  library(oligo)
  affy.raw <- read.celfiles(fullpaths)
  data.norm <- rma(affy.raw)
} else if(ALGORITHM.U == "MAS5" & !IS.OLIGO) {
  library(affy)
  if(COMPRESSED)
    stop("not implemented")

  data.all <- c()
  for(i1 in seq(1, length(fullpaths), MAS5.BATCH.SIZE)) {
    i2 <- min(i1+MAS5.BATCH.SIZE-1, length(fullpaths))
    print(c(i1, i2))
    x <- ReadAffy(filenames=(fullpaths[i1:i2]), compress=COMPRESSED)
    x.mas5 <- mas5(x, normalize=FALSE)
    data.all <- cbind(data.all, exprs(x.mas5))
  }
  data.norm <- new("ExpressionSet", exprs=data.all)
  data.norm <- affy.scalevalue.exprSet(data.norm)
} else if(ALGORITHM.U == "MAS5" & IS.OLIGO) {
  stop("Not implemented")
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
PSID.H <- "Probe.Set.ID"
DESCRIPTION.H <- "Target.Description"
LLID.H <- "LocusLink"
SYMBOL.H <- "Gene.Symbol"

if(!any(names(affyannot) == PSID.H))
  PSID.H <- "probeset_id"
if(!any(names(affyannot) == LLID.H))
  LLID.H <- "Entrez.Gene"
if(!any(names(affyannot) == LLID.H))
  LLID.H <- "geneid"

if(!any(names(affyannot) == PSID.H))
  stop("Cannot find probe set ID")
if(!any(names(affyannot) == DESCRIPTION.H))
  stop("Cannot find description")
if(!any(names(affyannot) == LLID.H))
  stop("Cannot find Entrez Gene ID")
if(!any(names(affyannot) == SYMBOL.H))
  stop("Cannot find gene symbol")

# Allow mismatches in control probes.
x <- rep(FALSE, length(probeset.id))
x[grep("^ERCC-", probeset.id)] <- TRUE
x[grep("^AFFX-BkGr-", probeset.id)] <- TRUE
I.control <- x

I <- match(probeset.id, affyannot[[PSID.H]])
if(any(is.na(I) & !I.control)) stop("mismatch")
affyannot <- affyannot[I,]
if(any(probeset.id[!I.control] != affyannot[[PSID.H]][!I.control])) stop("bug")

description <- affyannot[[DESCRIPTION.H]]
locuslink.id <- affyannot[[LLID.H]]
gene.symbol <- affyannot[[SYMBOL.H]]

description[is.na(description)] <- ""
locuslink.id[is.na(locuslink.id)] <- ""
gene.symbol[is.na(gene.symbol)] <- ""

# Write the annotated data.
data.out <- cbind(probeset.id, description, locuslink.id, gene.symbol, M)
colnames(data.out) <- c(
  "Probe.Set.ID", "Description", "LocusLink", "Gene.Symbol", sample.names)
my.write(data.out, FILE.NORM, col.names=TRUE)


# Clean up the TMP file.
unlink(FILE.TMP)
