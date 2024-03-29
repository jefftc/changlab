# clean.genes
#
# write.gmt
# read.gmt
# write.gmx

# A geneset is a list with members:
# NAME
# DESCRIPTION
# GENES


.format.gmt <- function(genesets) {
  maxlen <- max(unlist(lapply(genesets, function(x) length(x$GENES))))
  data.out <- c()
  for(gs in genesets) {
    description <- gs$DESCRIPTION
    if(description == "")
      description <- "na"
    x <- c(gs$GENES, rep("", maxlen-length(gs$GENES)))
    x <- c(gs$NAME, description, x)
    data.out <- rbind(data.out, x)
  }
  data.out
}

clean.genes <- function(genes) {
  x <- genes
  x <- x[x != "---"]
  x <- x[x != ""]
  x <- x[!is.na(x)]

  # Sort numerically, of possible.
  ow <- options("warn")
  options(warn=-1)
  nx <- as.numeric(x)
  options(ow)
  if(!any(is.na(nx)))
    x <- nx
  x <- sort(x)
  x <- x[!duplicated(x)]
  as.character(x)
}

write.gmt <- function(filename, genesets) {
  data.out <- .format.gmt(genesets)
  write.table(data.out, file=filename, quote=FALSE, sep="\t",
    row.names=FALSE, col.names=FALSE)
}

read.gmt <- function(filename) {
  if(length(grep("\\.gz$", filename, perl=TRUE, ignore.case=TRUE)))
    filename <- gzfile(filename)
  data <- read.delim(filename, header=FALSE, as.is=TRUE, comment.char="", 
    quote="", skip=0, colClasses=NA)

  genesets <- list()
  for(i in 1:nrow(data)) {
    x <- as.character(data[i,])
    NAME <- x[1]
    DESCRIPTION <- x[2]
    GENES <- x[3:length(x)]
    GENES <- GENES[GENES != ""]
    x <- list(NAME=NAME, DESCRIPTION=DESCRIPTION, GENES=GENES)
    genesets[[NAME]] <- x
  }
  genesets
}

write.gmx <- function(filename, genesets) {
  data.out <- .format.gmt(genesets)
  data.out <- t(data.out)
  write.table(data.out, filename, quote=FALSE, sep="\t",
    row.names=FALSE, col.names=FALSE)
}
