### Calculate adaptive confidence interval limits
### Functions according to definitions in Wassmer and Brannath

## Define generic combination function for Inverse Normal and MCT (2-stage)
comb <- function(x, y, weight) {
  lw <- length(weight)
  1 - pnorm(pmax.int(
    sqrt(weight[1]) * qnorm(1 - x) + sqrt(1 - weight[1]) * qnorm(1 - y),
    sqrt(weight[lw]) * qnorm(1 - x) + sqrt(1 - weight[lw]) * qnorm(1 - y)
  ))
}

## Define p-value as function of d (delta)
pd <- function(d, diff, sem, df, lt = FALSE) {
  pt((diff - d) / sem, df = df, lower.tail = lt)
}

## Define indicator function for equation (8.8)
ind <- function(x, y, d, weight, diff1, diff2, sem1, sem2, df1, df2, lt = FALSE) {
  comb(x, y, weight) <= comb(pd(d, diff1, sem1, df1, lt), 
                             pd(d, diff2, sem2, df2, lt), weight)
}

## Vectorized version of ind()
indv <- function(x, d, weight, diff1, diff2, sem1, sem2, df1, df2, lt = FALSE) {
  matrix(apply(x, 2, 
               function(z) comb(z[1], z[2], weight) <= 
                 comb(pd(d, diff1, sem1, df1, lt), pd(d, diff2, sem2, df2, lt), 
                      weight)), 
         ncol = ncol(x))
}

## Define function Qd from (8.8) as function of d (and subtract a)
Qd_a <- function(d, a1, a0, a, weight, diff1, diff2, sem1, sem2, df1, df2, 
                 lt = FALSE) {
  a1 + cubature::hcubature(f = indv, lowerLimit = c(a1, 0), 
                           upperLimit = c(a0, 1), d = d, weight = weight,
                           diff1 = diff1, diff2 = diff2, sem1 = sem1, sem2 = sem2,
                           df1 = df1, df2 = df2, lt = lt, maxEval = 10000, 
                           vectorInterface = TRUE)$integral - a
}

## Define function to calculate adaptive (1 - a) confidence region
## according to Section 8.2.1
adaptive_ci_limit <- function(diff1, diff2, sem1, sem2, df1, df2, 
                              a1, a0, a, weight, lower_bnd = TRUE) {
  # lower_bnd = TRUE creates lower bound l from one-sided interval (l, Inf)
  # lower_bnd = FALSE creates upper bound u from one-sided interval (-Inf, u)
  
  # Need root of Qd_a
  diff_tmp <- (df1 * diff1 + df2 * diff2) / (df1 + df2)
  search_int <- diff_tmp + c(-6, 6) * max(sem1, sem2)
  lt <- if (lower_bnd) FALSE else TRUE
  uniroot(f = Qd_a, interval = search_int, a1 = a1, a0 = a0, a = a, 
          weight = weight, diff1 = diff1, diff2 = diff2, sem1 = sem1,
          sem2 = sem2, df1 = df1, df2 = df2, lt = lt, tol = 1e-06)$root
}

## Define Median Unbiased Point Estimate according to Section 8.3.3
median_unbiased_pe <- function(diff1, diff2, sem1, sem2, df1, df2,
                               a1, a0, weight, lower_bnd = TRUE) {
  # Default setting via lower bound l from one-sided 50% interval (l, Inf)
  adaptive_ci_limit(diff1, diff2, sem1, sem2, df1, df2, a1, a0, 0.5, weight,
                    lower_bnd)
}

## Repeated confidence bounds for combination tests according to Section 8.2.2
repeated_ci <- function(diff1, diff2 = NULL, sem1, sem2 = NULL, df1, df2 = NULL,
                        a1, a2 = NULL, weight = NULL, stage = 1) {
  tol <- 1e-06
  if (stage == 1) {
    search_int <- diff1 + c(-6, 6) * sem1
    l <- uniroot(f = function(d) pd(d, diff1, sem1, df1, lt = FALSE) - a1,
                 interval = search_int, tol = tol)$root
    u <- uniroot(f = function(d) pd(d, diff1, sem1, df1, lt = TRUE) - a1,
                 interval = search_int, tol = tol)$root
  } else if (stage == 2) {
    stopifnot(!is.null(diff2), !is.null(sem2), !is.null(df2), !is.null(a2),
              !is.null(weight))
    diff_tmp <- (df1 * diff1 + df2 * diff2) / (df1 + df2)
    search_int <- diff_tmp + c(-6, 6) * max(sem1, sem2)
    l <- uniroot(f = function(d) comb(pd(d, diff1, sem1, df1, lt = FALSE), 
                                      pd(d, diff2, sem2, df2, lt = FALSE), 
                                      weight) - a2, 
                 interval = search_int, tol = tol)$root
    u <- uniroot(f = function(d) comb(pd(d, diff1, sem1, df1, lt = TRUE), 
                                      pd(d, diff2, sem2, df2, lt = TRUE), 
                                      weight) - a2, 
                 interval = search_int, tol = tol)$root
  } else {
    stop("Only stage equal to 1 or 2 implemented.")
  }
  if (l <= u) list(lower = l, upper = u) else NA
}