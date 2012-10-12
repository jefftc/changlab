# calc.km
# calc.km2               Different algorithm.
# calc.km.multi          Can do multiple groups.
# 
# plot.km
# plot.km.multi
# 
# write.km.prism
# write.km.multi.prism   Write out for multiple groups.

.calc.surv50 <- function(survival, dead) {
  library("survival")
  sf <- survfit(Surv(survival, dead) ~ rep("x", length(survival)))
  s <- summary(sf)
  surv.90 <- NA; surv.50 <- NA
  if(min(abs(s$surv-0.90)) < 0.05) {  # within 5%
    i <- which.min(abs(s$surv-0.90))
    surv.90 <- s$time[i]
  }
  if(min(abs(s$surv-0.50)) < 0.05) {  # within 5%
    i <- which.min(abs(s$surv-0.50))
    surv.50 <- s$time[i]
  }
  list(surv.90=surv.90, surv.50=surv.50)
}

.calc.km.curve <- function(survival, dead) {
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


# Split a vector of numeric values into groups.  The groups are split
# according to breakpoints, which is a vector of numbers from 0 to 1.
# For example, breakpoints of c(0.25, 0.50, 0.75) will split the
# values into 4 groups, the ones in the lower 25%, 25-50%, 50-75%, and
# greater than 75%.  breakpoints c(0.10) will split into two groups,
# 0-10%, and 10-100%.
group.by.value <- function(values, breakpoints) {
  # Make sure breakpoints are sorted, all >= 0 and <= 1.
  cutoffs <- sort(breakpoints)
  if((min(cutoffs) < 0) | (max(cutoffs > 1))) stop("bad breakpoints")
  values.r <- rank(values)
  ranks <- sapply(cutoffs, function(x) round(x*(length(values)-1))+1)
  group <- sapply(values.r, function(x) sum(x >= ranks))
  group
}

# group should be a vector that indicates which group each sample
# belongs to, e.g. c("LO", "LO", "MED", "MED", "HIGH") or c(0, 0, 1,
# 1, 2).
calc.km.multi <- function(survival, dead, group) {
  library("survival")
  if(length(survival) != length(dead)) stop("unaligned")
  if(length(survival) != length(group)) stop("unaligned")

  status <- factor(group)
  res <- coxph(Surv(survival, dead) ~ status, method="breslow")
  # rho=0 does log-rank test
  sd <- survdiff(Surv(survival, dead) ~ status, rho=0)
  df <- length(unique(group))-1
  p.value <- 1 - pchisq(res$score, df)
  hr <- exp(res$coefficients)

  surv <- list()
  for(g in unique(group)) {
    km <- .calc.surv50(survival[g==group], dead[g==group])
    surv[[as.character(g)]] <- km
  }

  list(p.value=p.value, hr=hr, surv=surv)
}

plot.km <- function(survival1, dead1, survival2, dead2, col1=NA, col2=NA, 
  main=NA, sub=NA, xlab="", ylab="") {
  if(is.na(col1))
    col1 <- "#000000"
  if(is.na(col2))
    col2 <- "#000000"

  if(length(survival1) != length(dead1)) stop("unaligned")
  if(length(survival2) != length(dead2)) stop("unaligned")
  survival <- c(survival1, survival2)
  dead <- c(dead1, dead2)
  xlim <- c(0, max(survival))
  ylim <- c(0, 1.00)

  km.1 <- .calc.km.curve(survival1, dead1)
  km.2 <- .calc.km.curve(survival2, dead2)

  plot(NA, type="n", axes=TRUE, xlim=xlim, ylim=ylim, 
    main=main, sub=sub, xlab=xlab, ylab=ylab)
  lines(km.1$surv.x, km.1$surv.y, col=col1)
  lines(km.2$surv.x, km.2$surv.y, col=col2)
  # Draw the censor lines.
  points(km.1$cens.x, km.1$cens.y, pch=15, cex=0.4)
  points(km.2$cens.x, km.2$cens.y, pch=15, cex=0.4)
}

plot.km.multi <- function(survival, dead, group, col=NA) {
  if(is.na(col))
    col <- list()
  if(length(survival) != length(dead)) stop("unaligned")
  if(length(survival) != length(group)) stop("unaligned")

  xlim <- c(0, max(survival))
  ylim <- c(0, 1.00)

  plot(NA, type="n", axes=TRUE, xlim=xlim, ylim=ylim, xlab="", ylab="")
  for(g in group) {
    co <- col[[as.character(g)]]
    if(is.null(co))
      co <- "#000000"
    km <- .calc.km.curve(survival[g==group], dead[g==group])
    lines(km$surv.x, km$surv.y, col=co)
    # Draw the censor lines.
    points(km$cens.x, km$cens.y, pch=15, cex=0.4)
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
write.km.prism.multi <- function(filename, survival, dead, group) {
  if(length(survival) != length(dead)) stop("unaligned")
  if(length(survival) != length(group)) stop("unaligned")

  all.groups <- sort(unique(group))
  data.out <- survival
  for(g in all.groups) {
    x <- rep("", length(dead))
    x[group==g] <- dead[group==g]
    data.out <- cbind(data.out, x)
  }
  data.out <- cbind(data.out, dead)
  colnames(data.out) <- c("Survival", all.groups, "All")
  write.table(data.out, filename, quote=FALSE, sep="\t",
    row.names=FALSE, col.names=TRUE)
}

