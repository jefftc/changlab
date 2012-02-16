load.and.merge.datasets <- function(file1, file2, map.file) {
  data.1 <- read.delim(file1,
    header=TRUE, as.is=TRUE, comment.char="", quote="")
  data.2 <- read.delim(file2,
      header=TRUE, as.is=TRUE, comment.char="", quote="")
  data.map <- read.delim(map.file, 
    header=FALSE, as.is=TRUE, comment.char="", quote="")
  data.1 <- data.1[data.map[,1],]
  data.2 <- data.2[data.map[,2],]
  list(data.1, data.2)
}
