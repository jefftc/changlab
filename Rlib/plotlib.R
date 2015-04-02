# green.shade
# red.shade
# rg.array.colors       Red/green array.
# by.array.colors       Blue/yellow array.
# rgb.colors            From red, green, blue colorwheel.
# ryb.colors            From red, yellow, blue colorwheel.
# matlab.colors         Default Matlab "jet" colors.
# matlab.hot.colors
# broad.colors
# bild.colors           Andrea's version of the "jet" colors.
# yahoo.weather.colors  Like Yahoo Weather Map.
# genespring.colors
# 
# my.lineplot
# my.lineplot2
# my.heatmap
# my.colorbar
# sample.colorbars      Print out all different color schemes.
# my.groupplot
# my.pcaplot
# my.boxplot
#
# INTERNAL FUNCTIONS
# matrix2color
# matrix2rgb


matrix2color <- function(cmatrix, pos) {
  # pos is [0, 1].  Returns r, g, b where each one is from [0, 1].
  if(is.nan(pos))  # this can happen if someone calls matlab.colors(1)
    pos <- 0.5
  breaks <- cmatrix[,1]
  i1 <- sum(pos >= breaks)
  x <- cmatrix[i1,2:4]
  if(i1 < nrow(cmatrix)) {
    i2 <- i1 + 1
    delta <- (pos - cmatrix[i1,1]) / (cmatrix[i2,1]-cmatrix[i1,1])
    x <- cmatrix[i1,2:4] + delta*(cmatrix[i2,2:4]-cmatrix[i1,2:4])
  }
  x/255
}

matrix2rgb <- function(cmatrix, pos) {
  # pos is [0, 1].  Returns color, e.g. "#003300".
  x <- matrix2color(cmatrix, pos)
  rgb(x[1], x[2], x[3], maxColorValue=1)
}

green.shade <- function(n) {
  if(n <= 0) stop("n must be greater than 0")
  cmatrix <- t(matrix(c(
    c(0.0,   0,   0, 0),
    c(1.0,   0, 255, 0)),
    4, 2))
  sapply((0:(n-1))/(n-1), function(x) matrix2rgb(cmatrix, x))
}

red.shade <- function(n) {
  if(n <= 0) stop("n must be greater than 0")
  cmatrix <- t(matrix(c(
    c(0.0,   0, 0, 0),
    c(1.0, 255, 0, 0)),
    4, 2))
  sapply((0:(n-1))/(n-1), function(x) matrix2rgb(cmatrix, x))
}

rg.array.colors <- function(n) {
  if(n <= 0) stop("n must be greater than 0")
  cmatrix <- t(matrix(c(
    c(0.0,   0, 255, 0),
    c(0.5,   0,   0, 0),
    c(1.0, 255,   0, 0)),
    4, 3))
  sapply((0:(n-1))/(n-1), function(x) matrix2rgb(cmatrix, x))
}

by.array.colors <- function(n) {
  if(n <= 0) stop("n must be greater than 0")
  cmatrix <- t(matrix(c(
    c(0.0,   0,   0, 255),
    c(0.5,   0,   0,   0),
    c(1.0, 255, 255,   0)),
    4, 3))
  sapply((0:(n-1))/(n-1), function(x) matrix2rgb(cmatrix, x))
}

rgb.colors <- function(n) {
  if(n <= 0) stop("n must be greater than 0")
  cmatrix <- t(matrix(c(
    c(0.00,   0,   0, 255),
    c(0.25,   0, 255, 255),
    c(0.50,   0, 255,   0),
    c(0.75, 255, 255,   0),
    c(1.00, 255,   0,   0)),
    4, 5))
  sapply((0:(n-1))/(n-1), function(x) matrix2rgb(cmatrix, x))
}

ryb.colors <- function(n) {
  if(n <= 0) stop("n must be greater than 0")
  cmatrix <- t(matrix(c(
    c(0.00,   0,   0, 255),
    c(0.25,   0, 255,   0),
    c(0.50, 255, 255,   0),
    c(0.75, 255, 128,   0),
    c(1.00, 255,   0,   0)),
    4, 5))
  sapply((0:(n-1))/(n-1), function(x) matrix2rgb(cmatrix, x))
}

matlab.colors <- function(n) {
  if(n <= 0) stop("n must be greater than 0")
  cmatrix <- t(matrix(c(
    c(0.000,   0,   0, 143),
    c(0.125,   0,   0, 255),
    c(0.250,   0, 127, 255),
    c(0.375,   0, 255, 255),
    c(0.500, 127, 255, 127),
    c(0.625, 255, 255,   0),
    c(0.750, 255, 127,   0),
    c(0.875, 255,   0,   0),
    c(1.000, 127,   0,   0)),
    4, 9))
  sapply((0:(n-1))/(n-1), function(x) matrix2rgb(cmatrix, x))
}

matlab.hot.colors <- function(n) {
  if(n <= 0) stop("n must be greater than 0")
  cmatrix <- t(matrix(c(
    c(0.000,  10,   0,   0),
    c(0.361, 255,   0,   0),
    c(0.377, 255,  10,   0),
    c(0.738, 255, 255,   0),
    c(1.000, 255, 255, 255)),
    4, 5))
  sapply((0:(n-1))/(n-1), function(x) matrix2rgb(cmatrix, x))
}

broad.colors <- function(n) {
  if(n <= 0) stop("n must be greater than 0")
  cmatrix <- t(matrix(c(
    c(0.000,  69,   0, 173),
    c(0.091,  39,   0, 209),
    c(0.182, 107,  88, 239),
    c(0.273, 136, 136, 255),
    c(0.364, 199, 193, 255),
    c(0.455, 213, 213, 255),
    c(0.545, 255, 192, 229),
    c(0.636, 255, 137, 137),
    c(0.727, 255, 112, 128),
    c(0.818, 255,  90,  90),
    c(0.909, 239,  64,  64),
    c(1.000, 214,  12,   0)),
    nrow=4))
  sapply((0:(n-1))/(n-1), function(x) matrix2rgb(cmatrix, x))
}

bild.colors <- function(n) {
  if(n <= 0) stop("n must be greater than 0")
  cmatrix <- t(matrix(c(
    c(0.000,  49,  50, 114),
    c(0.050,  61,  69, 137),
    c(0.100,  62,  84, 154),
    c(0.150,  67,  89, 160),
    c(0.200,  85, 108, 176),
    c(0.250, 115, 145, 201),
    c(0.300, 160, 205, 240),
    c(0.350, 180, 220, 243),
    c(0.400, 169, 216, 211),
    c(0.450, 160, 208, 164),
    c(0.500, 179, 213, 112),
    c(0.550, 203, 220,  61),
    c(0.600, 232, 231,  61),
    c(0.650, 255, 234,  47),
    c(0.700, 250, 180,  50),
    c(0.750, 243, 136,  54),
    c(0.800, 231,  80,  61),
    c(0.850, 218,  54,  55),
    c(0.900, 204,  55,  59),
    c(0.950, 160,  52,  52),
    c(1.000, 114,  39,  44)),
    4, 21))
  sapply((0:(n-1))/(n-1), function(x) matrix2rgb(cmatrix, x))
}

yahoo.weather.colors <- function(n) {
  if(n <= 0) stop("n must be greater than 0")
  cmatrix <- t(matrix(c(
    c(0.0, 255, 255, 255),
    c(0.1, 204, 255, 255),
    c(0.2, 153, 255, 255),
    c(0.3, 102, 204, 255),
    c(0.4,  84, 169, 255),
    c(0.5, 204, 255, 103),
    c(0.6, 255, 255, 103),
    c(0.7, 255, 204, 102),
    c(0.8, 255, 153, 102),
    c(0.9, 204, 102, 102),
    c(1.0, 209,  73,  73)),
    4, 11))
  sapply((0:(n-1))/(n-1), function(x) matrix2rgb(cmatrix, x))
}

genespring.colors <- function(n) {
  if(n <= 0) stop("n must be greater than 0")
  cmatrix <- t(matrix(c(
    c(0.0,   0,   0, 255),
    c(0.5, 255, 255,   0),
    c(1.0, 255,   0,   0)),
    4, 3))
  sapply((0:(n-1))/(n-1), function(x) matrix2rgb(cmatrix, x))
}



my.lineplot <- function(coords, xlim=NA, ylim=NA, col=NA, lwd=2) {
  # coords is a matrix where the columns are X_1, Y_1, ... X_N, Y_N.
  X <- coords[,seq(1, ncol(coords), 2)]
  Y <- coords[,seq(2, ncol(coords), 2)]
  if(ncol(X) != ncol(Y)) stop("X and Y unaligned")
  if(is.na(xlim))
    xlim <- c(min(X), max(X))
  if(is.na(ylim))
    ylim <- c(min(Y), max(Y))
  if(is.na(col))
    col <- matlab.colors(ncol(X))
  if(length(col) != ncol(X)) stop("colors unaligned")

  plot(NA, type="n", axes=TRUE, xlim=xlim, ylim=ylim, xlab="", ylab="")
  for(i in 1:ncol(X)) {
    x <- X[,i]; y <- Y[,i]
    lines(x, y, lwd=lwd, col=col[i])
    points(x, y, col=col[i], pch=".", cex=8)
  }
}

my.lineplot2 <- function(X, Y=NA, xlab=NA, ylab=NA, 
  col=NA, SPACING=NA, lwd=1, main=NA) {
  # X should be gene x samples matrix.
  # Y should be a vector of the line number for each row of X (1-based).
  # col is a vector of the colors for each gene.
  # By default, prints from the bottom (row 1) up.
  if((length(Y) == 1) && is.na(Y))
    Y <- 1:nrow(X)
  num.lines <- max(Y)
  if((length(col)==1) && is.na(col))
    col <- rep("#000000", nrow(X))
  if((length(ylab)==1) && is.na(ylab))
    ylab <- NA
  if((length(xlab)==1) && is.na(xlab))
    xlab <- NA
  if((length(SPACING) == 1) && is.na(SPACING))
    SPACING <- max(abs(X))

  xlim <- c(1, ncol(X))  
  ylim <- c(-SPACING, num.lines*SPACING)
  plot(NA, type="n", axes=FALSE, xlim=xlim, ylim=ylim, xlab="", ylab="", 
    main=main)
  for(i in 1:nrow(X)) {
    offset <- (Y[i]-1) * SPACING
    x <- as.numeric(X[i,]) + offset
    lines(x, lwd=lwd, col=col[i])
  }
  axis(1, at=1:ncol(X), labels=xlab)
  axis(2, at=seq(0, (num.lines-1)*SPACING, SPACING), labels=ylab)
  box()
}

normalize.one.mv <- function(x, M=0, V=1) {
  # Normalize a list of numbers so the mean is M and variance is V.
  M.0 <- mean(x)
  V.0 <- var(x)
  if(is.null(M))
    M <- M.0
  if(is.null(V))
    V <- V.0
  if(V.0 == 0)
    return(x-M.0+M)
  (x-M.0)*sqrt(V/V.0) + M
}

normalize.mv <- function(X, M=0, V=1) {
  t(apply(X, 1, function(x) normalize.one.mv(x, M, V)))
}

# matrix should contain values from 0 to 1.
my.heatmap <- function(matrix, col=rg.array.colors, xlab="", ylab="", 
  normalize=FALSE, scale=FALSE, cluster=FALSE) {
  # If normalize is TRUE, then will normalize to N(0, 1).  It can also
  # be a vector of (mean, variance) to specify the normalization
  # parameters.
  # If scale is TRUE, then will everything onto a 0 to 1 scale.  It
  # can also be a vector of (min, max, median) that indicates the
  # minimum, maximum, and median values that should correspond to 0,
  # 1, and 0.5.  If median is NA, then will not attempt to center the
  # median.
  # If cluster is TRUE, will cluster the rows and columns.
  if((length(normalize) == 1 && normalize == TRUE) || 
    (length(normalize) == 2)) {
    M <- 0; V <- 1
    if(length(normalize) == 2) {
      M <- normalize[1]; V <- normalize[2] 
    }
    m <- apply(matrix, 1, mean) - M
    matrix <- sweep(matrix, 1, m)
    matrix <- t(apply(matrix, 1, function(x) {
      V.0 <- var(x); M.0 <- mean(x); (x-M.0)*sqrt(V/V.0) + M.0 }))
  }

  if((length(scale) == 1 && scale == TRUE) || (length(scale) == 3)) {
    scale.min <- NA; scale.max <- NA; scale.med <- NA
    if(length(scale) != 1) {
      scale.min <- scale[1]; scale.max <- scale[2]; scale.med <- scale[3]
    }
    if(is.na(scale.min)) scale.min <- min(matrix)
    if(is.na(scale.max)) scale.max <- max(matrix)
    if(scale.max <= scale.min) stop("invalid scale parameters")
    matrix[matrix < scale.min] <- scale.min
    matrix[matrix > scale.max] <- scale.max
    if(!is.na(scale.med)) {
      # Center around 0, then scale so it's from -0.5 to +0.5.
      matrix <- matrix - scale.med
      x <- max(abs(c(min(matrix), max(matrix))))
      matrix <- matrix / x / 2 + 0.5
    } else {
      matrix <- matrix - min(matrix)
      matrix <- matrix / max(matrix)
    }
  }

  # col should be dark to bright.
  if(mode(col) == "character") {
    num.colors <- length(col)
  } else if(mode(col) == "list") {
    num.colors <- 256
    col <- rg.array.colors(num.colors)
  } else {
    num.colors <- 256
    col <- col(num.colors)
  }
  x.size <- 1
  y.size <- 1
  matrix <- as.matrix(matrix)

  # image treats the rows as x and columns as y.  So I need to
  # "rotate" the matrix 90 degrees clockwise.
  matrix <- matrix(t(matrix)[,nrow(matrix):1], ncol(matrix), nrow(matrix))

  # Weird.  If X11 is not running, then have to do it this way.
  # image puts row 1 on the bottom and column 1 on the right.  Flip
  # both the rows and columns.
  #matrix <- matrix[nrow(matrix):1,ncol(matrix):1]

  if(cluster) {
    h.r <- hclust(dist(matrix))
    h.c <- hclust(dist(t(matrix)))
    matrix <- matrix[h.r$order, h.c$order]
  }

  x <- x.size * 1:nrow(matrix)
  y <- y.size * 1:ncol(matrix)
  breaks <- (0:num.colors)/num.colors
  image(x, y, matrix, xlab=xlab, ylab=ylab, axes=FALSE, col=col, breaks=breaks)
  matrix
}

my.colorbar <- function(n=65, col=rg.array.colors) {
  # By default, Matlab plots n=65.
  I <- (n-1):0 / (n-1)
  M <- matrix(I, length(I), 1)
  my.heatmap(M, col=col)
}

.groupplot.calc.X <- function(
  Y, col.width=1.0, glyph.width=0.2, glyph.height=0.2, offset=0) {
  # x and y are the centers of the glyph.

  # Group each of the Y coordinates into bins of glyph_height.
  bin.height <- glyph.height
  bins <- sapply(Y, function(y) floor(y/bin.height))

  max.glyphs <- 1 + 2*col.width/glyph.width
  X.new <- c(); Y.new <- c(); I.new <- c()
  for(bin in unique(bins)) {
    I <- which(bins==bin)
    I <- I[order(Y[I])]
    Y.bin <- Y[I]

    delta <- glyph.width
    if(length(Y.bin) > max.glyphs)
      delta <- col.width / ((length(Y.bin)-1)/2)

    X.bin <- sapply(1:length(Y.bin), function(i) (i-1)/2*delta)
    if(length(X.bin) >= 2)
      X.bin[seq(2, length(X.bin), 2)] <- -X.bin[seq(2, length(X.bin), 2)]

    X.new <- c(X.new, X.bin)
    Y.new <- c(Y.new, Y.bin)
    I.new <- c(I.new, I)
  }

  X.new <- X.new + offset
  list(X.new, Y.new, I.new)
}

my.groupplot <- function(Ys, col.width=1.0, glyph.width=0.03, 
  glyph.height=0.03, col=NA) {
  # Ys is a list of vectors, where each vector is the Y-coordinates of
  # the points to plot.
  if(length(col.width) != 1) stop("invalid col.width")
  if(length(col)==1 && is.na(col))
    col <- lapply(Ys, function(y) rep("#000000", length(y)))

  num.groups <- length(Ys)

  offset = col.width
  usable.col.width = col.width/2.0 * 0.80

  Xs.new <- c(); Ys.new <- c(); col.new <- c()
  for(i in 1:num.groups) {
    Y <- Ys[[i]]
    x <- .groupplot.calc.X(
      Y, col.width=usable.col.width, glyph.width=glyph.width,
      glyph.height=glyph.height, offset=offset)
    X <- x[[1]]; Y <- x[[2]]; I <- x[[3]]
    Xs.new <- c(Xs.new, X)
    Ys.new <- c(Ys.new, Y)
    col.new <- c(col.new, col[[i]][I])
    offset <- offset + col.width
  }

  plot(Xs.new, Ys.new, pch=19, cex=1, frame.plot=FALSE, axes=FALSE, 
    col=col.new, xlab="", ylab="")
}

# Make a PCA plot of the samples (columns).
my.pcaplot <- function(X, K=3, d1=1, d2=2) {
  K <- max(c(K, d1, d2))
  S <- svd(X)
  Y <- t(S$u[,1:K]) %*% X
  plot(Y[d1,], Y[d2,], pch=".", cex=8, xlab="", ylab="")
  Y
}

sample.colorbars <- function() {
  M <- t(matrix(1:8, 4, 2))
  layout(M)
  col.fns <- c(rg.array.colors, by.array.colors, rgb.colors, ryb.colors, 
    matlab.colors, matlab.hot.colors, genespring.colors, yahoo.weather.colors)
  col.names <- c("Red/Green", "Yellow/Blue", "RGB", "RYB", 
    "Matlab Jet", "Matlab Hot", "GeneSpring", "Weather")
  # For some reason, I can't loop 1:length(col.fns).  col <-
  # col.fns[i] always gets me rg.array.colors.
  i <- 1
  for(col in col.fns) {
    my.colorbar(col=col)
    title(col.names[i])
    i <- i + 1
  }
}

my.boxplot <- function(X, labels=TRUE, col=NULL, main="", xlab="", ylab="", 
  leg=NULL, fill=NULL, cex.labels=1.25, cex.legend=1, vert.labels=TRUE,
  ylim=NULL, cex.lab=1.5, sub="", cex.sub=1.5, x.legend=NULL)
{
  las <- 0
  if(vert.labels)
    las <- 3
  at <- NULL
  if(length(labels) > 1)
    at <- 1:length(labels)
  if(is.null(x.legend))
    x.legend <- "bottomleft"

  lwd <- 2
  boxplot(X, col=col, xlab="", ylab="", axes=FALSE, pch=19, cex=1, ylim=ylim)
  box(lwd=lwd)
  axis(1, lwd=lwd, cex.axis=cex.labels, labels=labels, at=at, las=las)
  axis(2, lwd=lwd, cex.axis=1.5)
  title(main=main, xlab=xlab, ylab=ylab, cex.lab=cex.lab, cex.main=2.0, 
    sub=sub, cex.sub=cex.sub)
  if(length(leg))
    legend(x.legend, legend=leg, fill=fill, box.lwd=1.5, cex=cex.legend,
      bg="#FFFFFF", inset=0.05)
}
