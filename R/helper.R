critical.value.2stage <- function(alpha = 0.05, weight) {
  # TO DO: adapt to max.comb.test as well
  corr <- diag(1, nrow = 2, ncol = 2)
  corr[lower.tri(corr)] <- sqrt(weight)
  corr[upper.tri(corr)] <- sqrt(weight)
  
  # Could also use qmvnorm but is much slower
  f <- function(x) {
    pmvnorm(lower = rep.int(-Inf, nrow(corr)), upper = rep.int(x, nrow(corr)),
            mean = rep.int(0, nrow(corr)), corr = corr,
            algorithm = GenzBretz(maxpts = 100000, abseps = 1e-06)) - 
      (1 - alpha)
  }
  
  extint <- if (alpha > 0.5) "upX" else "no"
  z_crit_s1 <- uniroot(f, interval = c(0, qnorm(1 - alpha/2)),
                       extendInt = extint)$root
  list(cval = z_crit_s1, siglev = 1 - pnorm(z_crit_s1))
}