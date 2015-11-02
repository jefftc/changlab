# 151101  created

##################################################
# Extract parameters from the command line.
argv <- commandArgs()
if(length(argv) != 6)
  stop("Usage: sppscript.R --vanilla <treatment.bam> <control.bam> <fdr> <num.procs>")

TREATMENT.FILE <- argv[3]
CONTROL.FILE <- argv[4]
FDR <- argv[5]
NUM.PROCS <- argv[6]

# Check the parameters.
if(!file.exists(TREATMENT.FILE))
  stop("I cannot find treatment file.")
if(!file.exists(CONTROL.FILE))
  stop("I cannot find control file.")
if((FDR <= 0) | (FDR > 1))
  stop("Invalid FDR.")
if((NUM.PROCS < 1) | (NUM.PROCS > 100))
  stop("Invalid NUM.PROCS.")



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
# Run SPP.

library(spp)
library(snow)

cluster <- makeCluster(NUM.PROCS)

start <- proc.time()

chip.data <- read.bam.tags(TREATMENT.FILE)
input.data <- read.bam.tags(CONTROL.FILE)


## Choosing alignment quality, removing anomalies.

# Get binding info from cross-correlation profile.
# srange  possible range for the size of the protected region
#         Should be higher than tag length.  Making the upper boundary 
#         too high will increase calculation time.
# bin     bin tags within the specified number of basepairs to speed up
#         calculation.  increasing bin size decreases the accuracy of 
#         the determined parameters
binding.characteristics <- get.binding.characteristics(
  chip.data, srange=c(50, 500), bin=5, cluster=cluster)
print(paste("binding peak separation distance =", 
  binding.characteristics$peak$x))

# plot cross-correlation profile
pdf(file="crosscorrelation.pdf", width=5, height=5)
par(mar=c(3.5,3.5,1.0,0.5), mgp=c(2,0.65,0), cex=0.8)
plot(binding.characteristics$cross.correlation, type='l',
  xlab="strand shift", ylab="cross-correlation")
abline(v=binding.characteristics$peak$x, lty=2, col=2)
dev.off()

# select informative tags based on the binding characteristics
chip.data <- select.informative.tags(chip.data, binding.characteristics)
input.data <- select.informative.tags(input.data, binding.characteristics)

# restrict or remove singular positions with very high tag counts
chip.data <- remove.local.tag.anomalies(chip.data)
input.data <- remove.local.tag.anomalies(input.data)

## Calculating genome-wide tag density and tag enrichment/depletion
## profiles.

# output smoothed tag density (subtracting re-scaled input) into a
# WIG file note that the tags are shifted by half of the peak
# separation distance
# Use scale.by.dataset.size=T option to normalize the tag density by
# the total dataset size (to make it comparable across samples
# sequenced at different sequencing depth. note: this option is
# typically used without the background subtraction).
print("get.smoothed.tag.density")
tag.shift <- round(binding.characteristics$peak$x/2)
smoothed.density <- get.smoothed.tag.density(
  chip.data, control.tags=input.data, bandwidth=200, step=100, 
  tag.shift=tag.shift)
writewig(smoothed.density, "density.wig",
  "Smoothed, background-subtracted tag density")
rm(smoothed.density)

# To provide a rough estimate of the enrichment profile (i.e. ChIP
# signal over input).
print("smoothed.enrichment")
smoothed.enrichment.estimate <- get.smoothed.enrichment.mle(
  chip.data, input.data, bandwidth=200, step=100, tag.shift=tag.shift)
writewig(smoothed.enrichment.estimate, "enrichment.wig",
  "Smoothed maximum likelihood log2 enrichment estimate")

# output conservative enrichment estimates
# alpha specifies significance level at which confidence intervals
# will be estimated
print("estimate enrichment")
enrichment.estimates <- get.conservative.fold.enrichment.profile(
  chip.data, input.data, fws=500, step=100, alpha=0.01)
x1 <- "Conservative fold-enrichment/depletion estimates"
x2 <- "shown on log2 scale"
x <- sprintf("%s %s", x1, x2)
writewig(enrichment.estimates, "enrichment.estimates.wig", x)
rm(enrichment.estimates)

# Broad regions of enrichment for a specified scale can be quickly
# identified
print("find broad enrichment clusters")
broad.clusters <- get.broad.enrichment.clusters(
  chip.data, input.data, window.size=1e3, z.thr=3,
  tag.shift=round(binding.characteristics$peak$x/2))
# write out in broadPeak format
write.broadpeak.info(broad.clusters, "broadPeak")

# use WTD method to call binding positions, using FDR of 1% and a
# window size estimated by the binding.characteristics
# binding detection parameters
# desired FDR (1%). Alternatively, an E-value can be supplied to the
# method calls below instead of the fdr parameter
print("calling binding positions")
# the binding.characteristics contains the optimized half-size for
# binding detection window
detection.window.halfsize <- binding.characteristics$whs
  
# determine binding positions using wtd method
bp <- find.binding.positions(signal.data=chip.data, control.data=input.data,
  fdr=FDR, whs=detection.window.halfsize, cluster=cluster)
# Alternatively, the binding positions can be determined using MTC
#  method (referred here as lwcc)
#bp <- find.binding.positions(signal.data=chip.data, control.data=input.data,
#  fdr=fdr, method=tag.lwcc, whs=detection.window.halfsize, cluster=cluster)
x <- sum(unlist(lapply(bp$npl,function(d) length(d$x))))
print(paste(fdr, "detected", x, "peaks"))
print((proc.time()-start)[3])

# output detected binding positions
output.binding.results(bp, "binding.positions.txt")
  
# Broader regions of enrichment associated with the determined peaks
# can be determined using the following:
print("find narrow peaks")
bp <- add.broad.peak.regions(
  chip.data, input.data, bp, window.size=1000, z.thr=3)
# output using narrowPeak format
write.narrowpeak.binding(bp, "narrowPeak")

### determine MSER
### note: this will take approximately 10-15x the amount of time the
### initial binding detection did
### The saturation criteria here is 99% consistency in the set of
### binding positions when adding 1e5 tags.
### To ensure convergence the number of subsampled chains (n.chains)
### should be higher (80)
##print("determine MSER")
### Get error here:
### Error in if (me <= menr) { : missing value where TRUE/FALSE needed
##mser <- get.mser(chip.data, input.data, step.size=1e5, test.agreement=0.99,
##  n.chains=8, cluster=cluster, fdr=fdr, method=tag.wtd, 
##  whs=detection.window.halfsize)
##print(paste("MSER at a current depth is",mser))
##
### interpolate MSER dependency on tag count
### note: this requires considerably more calculations than the
### previous steps (~ 3x more than the first MSER calculation)
### Here we interpolate MSER dependency to determine a point at which
### MSER of 2 is reached
### The interpolation will be based on the difference in MSER at the
### current depth, and a depth at 5e5 fewer tags (n.steps=6);
### evaluation of the intermediate points is omitted here to speed up
### the calculation (excluded.steps parameter)
### A total of 7 chains is used here to speed up calculation, whereas
### a higher number of chains (50) would give good convergence
##msers <- get.mser.interpolation(chip.data, input.data, step.size=1e5, 
##  test.agreement=0.99, target.fold.enrichment=2, n.chains=7, n.steps=6,
##  cluster=cluster, fdr=fdr, method=tag.wtd, whs=detection.window.halfsize)
##
##x <- round(unlist(lapply(msers,function(x) x$prediction))/1e6, 5)
##print(paste("predicted sequencing depth =", x, " million tags"))
