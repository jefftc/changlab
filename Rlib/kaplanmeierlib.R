# calc.km
# calc.km2   Different algorithm.
# plot.km
# write.km.prism

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

  list(p.value=p.value)
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
  list(p.value=p.value, hr=hr)
}

plot.km <- function(survival1, dead1, survival2, dead2, col1=NA, col2=NA) {
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

  km.1 <- calc.km.curve(survival1, dead1)
  km.2 <- calc.km.curve(survival2, dead2)

  plot(NA, type="n", axes=TRUE, xlim=xlim, ylim=ylim, xlab="", ylab="")
  lines(km.1$surv.x, km.1$surv.y, col=col1)
  lines(km.2$surv.x, km.2$surv.y, col=col2)
  # Draw the censor lines.
  points(km.1$cens.x, km.1$cens.y, pch=15, cex=0.4)
  points(km.2$cens.x, km.2$cens.y, pch=15, cex=0.4)
}

write.km.prism <- function(filename, survival1, dead1, survival2, dead2) {
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
  write.table(data.out, filename, quote=FALSE, sep="\t",
    row.names=FALSE, col.names=FALSE)
}

calc.km.curve <- function(survival, dead) {
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
