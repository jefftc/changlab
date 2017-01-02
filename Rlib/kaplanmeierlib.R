# calc.km
# calc.km2               Different algorithm.
# calc.km.multi          Can do multiple groups.
# 
# plot.km
# plot.km.multi
# 
# write.km.prism
# write.km.prism.multi   Write out for multiple groups.
#
# group.by.value


# Likely to get warning if only 1 member per group.
# Maybe should set a specific error message in this case.
# Warning message:
# In fitter(X, Y, strats, offset, init, control, weights = weights,  :
#   Loglik converged before variable  1 ; beta may be infinite. 

# Estimate the survival time based on the closest survival.
.est.survival <- function(surv.table, survival) {
  # surv.table
  # time   survival
  # 25.2   0.75
  # 38.4   0.50
  # 40.8   0.25
  # 50.4   0.00

  time <- NA

  # Look for the closest one.  If one found within 5%, return it.
  if(min(abs(surv.table$surv-survival)) < 0.05) {  # within 5%
    i <- which.min(abs(surv.table$surv-survival))
    time <- surv.table$time[i]
  } 
  # If the survival if higher than the highest survival, then
  # extrapolate the time.
  # No.  Should use next lower time.
  #else if(survival > surv.table$surv[1]) {
  #  surv1 <- surv.table$surv[1]
  #  time1 <- surv.table$time[1]
  #  frac <- (1.0-survival)/(1.0-surv1)
  #  time <- frac * time1
  #} 
  # If the survival if lower than the lowest survival, then
  # the time is indefinite.
  else if(survival < surv.table$surv[length(surv.table$surv)]) {
    time <- NA
  }
  # Otherwise, the survival goes to the next lowest time.  E.g. if
  # looking for 60% survival, use the 50% survival time.
  else {
    I <- which(surv.table$surv < survival)
    if(length(I) < 1) stop("broken")
    time <- surv.table$time[I[1]]
  }
  
  time
}

.calc.surv50 <- function(survival, dead) {
  library("survival")
  surv.50 <- NA
  surv.90 <- NA
  sf <- survfit(Surv(survival, dead) ~ rep("x", length(survival)))
  s <- summary(sf)
  if(length(s$surv) > 0) {
    surv.90 <- .est.survival(s, 0.90)
    surv.50 <- .est.survival(s, 0.50)
  }
  list(surv.90=surv.90, surv.50=surv.50)
}

.calc.km.curve <- function(survival, dead) {
  if(any(is.na(survival))) stop("na in survival")
  if(any(is.na(dead))) stop("na in dead")
  O <- order(survival)
  survival <- survival[O]; dead <- dead[O]

  surv.x <- c(0)      # Start with 100% survival at time 0.
  surv.y <- c(1.00)
  cens.x <- c()
  cens.y <- c()

  for(i in 1:length(survival)) {
    if(dead[i]) {
      x <- survival[i]
      remaining.patients <- length(survival) - i + 1
      # Percent of the remaining patients have now died.
      perc <- 1/(remaining.patients)
      # Calculate percent survival.
      y <- surv.y[length(surv.y)] * (1-perc)
      surv.x <- c(surv.x, x, x)
      surv.y <- c(surv.y, surv.y[length(surv.y)], y)
    } else {
      # Censor this value.
      cens.x <- c(cens.x, survival[i])
      cens.y <- c(cens.y, surv.y[length(surv.y)])
    }
  }

  # If the last point was censored, then extend the survival line all
  # the way out to that patient.
  if(length(cens.x) && cens.x[length(cens.x)] > surv.x[length(surv.x)]) {
    x <- cens.x[length(cens.x)]
    y <- surv.y[length(surv.y)]
    surv.x <- c(surv.x, x)
    surv.y <- c(surv.y, y)
  }

  list(surv.x=surv.x, surv.y=surv.y, cens.x=cens.x, cens.y=cens.y)
}

##stratify.by.expression <- function(num.groups, expression, survival, dead) {
##  if(!length(expression)) 
##    stop("missing expression data")
##  if(length(expression) != length(survival))
##    stop("unaligned")
##  if(length(expression) != length(dead))
##    stop("unaligned")
##  if((num.groups < 1) | (num.groups > length(expression)))
##    stop("invalid num.groups")
##}

calc.km <- function(survival1, dead1, survival2, dead2) {
  # survival1 is a a vector of matrix times.  dead1 is 1/0 where 1
  # means the patient is dead and 0 otherwise.  Returns list of
  # p.value.
  library("survival")
  if(length(survival1) != length(dead1)) stop("unaligned")
  if(length(survival2) != length(dead2)) stop("unaligned")
  survival <- c(survival1, survival2)
  dead <- c(dead1, dead2)
  x <- c(rep("1", length(survival1)), rep("0", length(survival2)))
  status <- factor(x)
  sd <- survdiff(Surv(survival, dead) ~ status)
  p.value <- 1-pchisq(sd$chisq, 1)

  x1 <- .calc.surv50(survival1, dead1)
  x2 <- .calc.surv50(survival2, dead2)
  list(p.value=p.value, 
    surv1.50=x1$surv.50, surv1.90=x1$surv.90, 
    surv2.50=x2$surv.50, surv2.90=x2$surv.90)
}

calc.km2 <- function(survival1, dead1, survival2, dead2) {
  # survival1 is a a vector of matrix times.  dead1 is 1/0 where 1
  # means the patient is dead and 0 otherwise.  Returns list of
  # p.value.
  library("survival")
  if(length(survival1) != length(dead1)) stop("unaligned")
  if(length(survival2) != length(dead2)) stop("unaligned")
  survival <- c(survival1, survival2)
  dead <- c(dead1, dead2)
  x <- c(rep("1", length(survival1)), rep("0", length(survival2)))
  status <- factor(x)

  res <- coxph(Surv(survival, dead) ~ status, method="breslow")
  p.value <- 1 - pchisq(res$score, 1)
  hr <- exp(res$coefficients)
  #summary(res)

  x1 <- .calc.surv50(survival1, dead1)
  x2 <- .calc.surv50(survival2, dead2)
  list(p.value=p.value, hr=hr,
    surv1.50=x1$surv.50, surv1.90=x1$surv.90, 
    surv2.50=x2$surv.50, surv2.90=x2$surv.90)
}


# group should be a vector that indicates which group each sample
# belongs to, e.g. c("LO", "LO", "MED", "MED", "HIGH") or c(0, 0, 1,
# 1, 2).
calc.km.multi <- function(survival, dead, group) {
  library("survival")
  if(length(survival) != length(dead)) stop("unaligned")
  if(length(survival) != length(group)) stop("unaligned")

  group.all <- sort(unique(group[!is.na(group)]))

  status <- factor(group)
  res <- coxph(Surv(survival, dead) ~ status, method="breslow")
  # rho=0 does log-rank test
  sd <- survdiff(Surv(survival, dead) ~ status, rho=0)
  df <- length(group.all)-1
  p.value <- 1 - pchisq(res$score, df)
  hr <- exp(res$coefficients)

  num.samples <- list()
  for(g in group.all)
    num.samples[[as.character(g)]] <- sum(!is.na(group) & (group == g))

  surv <- list()
  for(g in group.all) {
    km <- .calc.surv50(survival[g==group], dead[g==group])
    surv[[as.character(g)]] <- km
  }

  list(p.value=p.value, hr=hr, surv=surv, num.samples=num.samples)
}

plot.km <- function(survival1, dead1, survival2, dead2, col1=NA, col2=NA, 
  main=NA, sub=NA, xlab="", ylab="", cex.main=NULL, name1=NA, name2=NA) {
  if(is.null(cex.main))
    cex.main <- 2.0
  if(is.na(col1))
    col1 <- "#000000"
  if(is.na(col2))
    col2 <- "#000000"

  if(length(survival1) != length(dead1)) stop("unaligned")
  if(length(survival2) != length(dead2)) stop("unaligned")
  survival <- c(survival1, survival2)
  dead <- c(dead1, dead2)
  xlim <- c(0, max(survival))
  #ylim <- c(0, 1.00)
  ylim <- c(0, 100)

  km.1 <- .calc.km.curve(survival1, dead1)
  km.2 <- .calc.km.curve(survival2, dead2)

  lwd <- 2
  main.line <- NA
  cex.sub <- 1.5
  cex.lab <- 2.5
  sub.line <- NA
  lwd.line <- 3
  cex.censor <- 0.8
  cex.legend <- 2.0

  plot(NA, type="n", axes=FALSE, xlim=xlim, ylim=ylim, xlab="", ylab="")
  usr <- par("usr")
  rect(usr[1], usr[3], usr[2], usr[4], col="#FFFFFF", border=NA)
  #box(lwd=lwd)
  axis(1, lwd=lwd, cex.axis=1.5)
  axis(2, lwd=lwd, cex.axis=1.5)
  title(main=main, cex.main=cex.main, line=main.line)
  title(xlab=xlab, ylab=ylab, cex.lab=cex.lab)
  title(sub=sub, cex.sub=cex.sub, col.sub="#A60400", line=sub.line)

  #plot(NA, type="n", axes=TRUE, xlim=xlim, ylim=ylim, 
  #  main=main, sub=sub, xlab=xlab, ylab=ylab)
  lines(km.1$surv.x, km.1$surv.y*100, col=col1, lwd=lwd.line)
  lines(km.2$surv.x, km.2$surv.y*100, col=col2, lwd=lwd.line)
  # Draw the censor lines.
  points(km.1$cens.x, km.1$cens.y*100, pch=15, cex=cex.censor)
  points(km.2$cens.x, km.2$cens.y*100, pch=15, cex=cex.censor)

  if(!is.na(name1) & !is.na(name2)) {
    leg <- c(name1, name2)
    fill <- c(col1, col2)
    legend("bottomright", legend=leg, fill=fill, bty="n", 
      box.lwd=1.5, cex=cex.legend, inset=0.05)
  }
}

# col is a list of NAME -> color (e.g. "#FF0000")
plot.km.multi <- function(survival, dead, group, col=NA, 
  main="", cex.main=NULL, main.line=NA, xlab="", ylab="", sub="", 
  cex.sub=NULL, sub.line=NA, cex.legend=NULL, x.legend=NULL) {
  if(is.null(cex.main))
    cex.main <- 2.0
  if(is.null(cex.sub))
    cex.sub <- 1.5
  if(is.null(cex.legend))
    cex.legend <- 1.5
  if(is.null(x.legend))
    x.legend <- "bottomleft"
  if(all(is.na(col)))
    col <- list()
  if(length(survival) != length(dead)) stop("unaligned")
  if(length(survival) != length(group)) stop("unaligned")

  # What is this for?
  #km <- calc.km.multi(survival, dead, group)

  xlim <- c(0, max(survival))
  ylim <- c(0, 100)

  lwd <- 2
  plot(NA, type="n", axes=FALSE, xlim=xlim, ylim=ylim, xlab="", ylab="")
  usr <- par("usr")
  rect(usr[1], usr[3], usr[2], usr[4], col="#FFFFFF")
  box(lwd=lwd)
  axis(1, lwd=lwd, cex.axis=1.5)
  axis(2, lwd=lwd, cex.axis=1.5)
  title(main=main, cex.main=cex.main, line=main.line)
  title(xlab=xlab, ylab=ylab, cex.lab=1.5)
  title(sub=sub, cex.sub=cex.sub, col.sub="#A60400", line=sub.line)

  all.groups <- sort(unique(group[!is.na(group)]))
  for(g in all.groups) {
    co <- col[[as.character(g)]]
    if(is.null(co))
      co <- "#000000"
    I <- !is.na(group) & (g == group)
    km <- .calc.km.curve(survival[I], dead[I])
    lines(km$surv.x, km$surv.y*100, col=co, lwd=lwd)
    # Draw the censor lines.
    points(km$cens.x, km$cens.y*100, pch=15, cex=0.8)
  }
  if(!all(is.na(col))) {
    leg <- names(col)
    fill <- sapply(leg, function(x) col[[x]])
    # Only label the lines if it exists on the plot.
    #I <- !is.na(match(leg, group))
    #leg <- leg[I]
    #fill <- fill[I]
    legend(x=x.legend, legend=leg, fill=fill, box.lwd=1.5, cex=cex.legend, 
      inset=0.05)
  }
}

write.km.prism <- function(filename, class1, survival1, dead1, 
  class2, survival2, dead2) {
  # Write out results for Prism.
  if(length(survival1) != length(dead1)) stop("unaligned")
  if(length(survival2) != length(dead2)) stop("unaligned")
  survival <- c(survival1, survival2)
  dead <- c(dead1, dead2)

  x.value <- survival
  a.value <- c(dead1, rep("", length(dead2)))
  b.value <- c(rep("", length(dead1)), dead2)
  c.value <- dead
  data.out <- cbind(x.value, a.value, b.value, c.value)
  colnames(data.out) <- c("Survival", class1, class2, "Both")
  write.table(data.out, filename, quote=FALSE, sep="\t",
    row.names=FALSE, col.names=TRUE)
}

# Write out results for Prism.
write.km.prism.multi <- function(filename, survival, dead, group, sample=NA) {
  if(length(survival) != length(dead)) stop("unaligned")
  if(length(survival) != length(group)) stop("unaligned")
  if(!all(is.na(sample))) {
    if(length(sample) != length(group)) stop("unaligned")
  }

  all.groups <- sort(unique(group[!is.na(group)]))
  data.out <- survival
  for(g in all.groups) {
    x <- rep("", length(dead))
    I <- !is.na(group) & (group == g)
    x[I] <- dead[I]
    data.out <- cbind(data.out, x)
  }
  data.out <- cbind(data.out, dead)
  colnames(data.out) <- c("Survival", all.groups, "All")
  if(!all(is.na(sample))) {
    data.out <- cbind(sample, data.out)
    colnames(data.out)[1] <- "Sample"
  }
  write.table(data.out, filename, quote=FALSE, sep="\t",
    row.names=FALSE, col.names=TRUE)
}

# Split a vector of numeric values into groups.  The groups are split
# according to breakpoints, which is a vector of numbers from 0 to 1.
# For example, breakpoints of c(0.25, 0.50, 0.75) will split the
# values into 4 groups, the ones in the lower 25%, 25-50%, 50-75%, and
# greater than 75%.  breakpoints c(0.10) will split into two groups,
# 0-10%, and 10-100%.
# 
# Groups are specified by numbers starting from 0.  Group 0 are the
# expression values lower than breakpoints[1].  So if breakpoints[1]
# is 0, then the lowest index will be 1.
group.by.value <- function(values, breakpoints) {
  # Make sure breakpoints are sorted, all >= 0 and <= 1.
  cutoffs <- sort(breakpoints)
  if((min(cutoffs) < 0) | (max(cutoffs > 1))) stop("bad breakpoints")
  values.r <- rank(values)
  ranks <- sapply(cutoffs, function(x) round(x*(length(values)-1))+1)
  group <- sapply(values.r, function(x) sum(x >= ranks))
  group
}
