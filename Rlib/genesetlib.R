# read.gmt
# write.gmt

# A geneset is a list with members:
# NAME
# DESCRIPTION
# GENES


write.gmt <- function(filename, genesets) {
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
  write.table(data.out, filename, quote=FALSE, sep="\t",
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
