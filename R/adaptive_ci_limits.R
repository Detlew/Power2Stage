### Calculate adaptive confidence interval limits

# Define generic combination function
comb <- function(x, y, weight) {
  lw <- length(weight)
  1 - pnorm(pmax.int(
    sqrt(weight[1]) * qnorm(1 - x) + sqrt(1 - weight[1]) * qnorm(1 - y),
    sqrt(weight[lw]) * qnorm(1 - x) + sqrt(1 - weight[lw]) * qnorm(1 - y)
  ))
  # What should be returned in case of x = 0 and y = 1 (or vice versa)?
}

# Define p-value as function of delta
pd <- function(d, diff, sem, df, lt = FALSE) {
  pt((diff - d) / sem, df = df, lower.tail = lt)
}

# Define vectorized indicator function for equation (8.8) in Wassmer + Brannath
ind <- function(x, y, d, weight, diff1, diff2, sem1, sem2, df1, df2, lt = FALSE) {
  comb(x, y, weight) <= comb(pd(d, diff1, sem1, df1, lt), 
                             pd(d, diff2, sem2, df2, lt), weight)
}

# vectorized version
indv <- function(x, d, weight, diff1, diff2, sem1, sem2, df1, df2, lt = FALSE) {
  # x generic variable
  # d = delta
  # lt = lower.tail
  matrix(apply(x, 2, 
               function(z) comb(z[1], z[2], weight) <= 
                 comb(pd(d, diff1, sem1, df1, lt), pd(d, diff2, sem2, df2, lt), 
                      weight)), 
         ncol = ncol(x))
}

# Define function Qd from (8.8) as function of d (Wassmer + Brannath)
Qd <- function(d, a1, a0, a, weight, diff1, diff2, sem1, sem2, df1, df2, 
               lt = FALSE) {
  # Note: cubature::hcubature seems to "run forever" if maxEval = 0 is set
  #       The reason does not seem to be a case where x = 0 and y = 0
  # Maybe use pracma::integral2 instead (with singular = TRUE argument?)
  a1 + cubature::hcubature(f = indv, lowerLimit = c(a1, 0), 
                           upperLimit = c(a0, 1), d = d, weight = weight,
                           diff1 = diff1, diff2 = diff2, sem1 = sem1, sem2 = sem2,
                           df1 = df1, df2 = df2, lt = lt, maxEval = 5000, 
                           vectorInterface = TRUE)$integral - a
  #a1 + pracma::integral2(fun = ind, xmin = a1, ymin = 0, xmax = a0, ymax = 1,
  #                       reltol = 1e-6, maxlist = 10000, singular = FALSE, 
  #                       d = d, lt = lt, weight = weight, diff1 = diff1, 
  #                       diff2 = diff2, sem1 = sem1, sem2 = sem2, df1 = df1, 
  #                       df2 = df2)$Q - a
}

adaptive_ci_limit <- function(diff1, diff2, sem1, sem2, df1, df2, 
                              a1, a0, a, weight, lower_bnd = TRUE) {
  # diff1 and diff2 are observed trt differences
  # n1 and n2 are actual sample sizes (not pre-defined) of stage 1 and 2, resp.
  # lower_bnd = TRUE creates lower bound l from one-sided interval (l, Inf)
  # lower_bnd = FALSE creates upper bound u from one-sided interval (-Inf, u)
  
  # Need root of Qd
  diff_tmp <- (df1 * diff1 + df2 * diff2) / (df1 + df2)
  search_int <- diff_tmp + c(-6, 6) * max(sem1, sem2)
  lt <- if (lower_bnd) FALSE else TRUE
  uniroot(f = Qd, interval = search_int, a1 = a1, a0 = a0, a = a, 
          weight = weight, diff1 = diff1, diff2 = diff2, sem1 = sem1,
          sem2 = sem2, df1 = df1, df2 = df2, lt = lt, tol = 1e-06)$root
}

# See section 8.3.3 in Brannath + Wassmer
median_unbiased_pe <- function(diff1, diff2, sem1, sem2, df1, df2,
                               a1, a0, weight, lower_bnd = TRUE) {
  # Default setting is to get median unbiased point estimate via lower bound l
  # from one-sided 50% interval (l, Inf)
  adaptive_ci_limit(diff1, diff2, sem1, sem2, df1, df2, a1, a0, 0.5, weight,
                    lower_bnd)
}

# Repeated confidence bounds for combination tests
# see section 8.2.2 in Brannath + Wassmer
repeated_ci <- function(diff1, diff2 = NULL, sem1, sem2 = NULL, df1, df2 = NULL,
                        a1, a2 = NULL, weight = NULL, stage = 1) {
  if (stage == 1) {
    search_int <- diff1 + c(-6, 6) * sem1
    # lower bound l from one-sided interval (l, Inf)
    l <- uniroot(f = function(d) pd(d, diff1, sem1, df1, lt = FALSE) - a1,
                 interval = search_int, tol = 1e-06)$root
    # upper bound u from one-sided interval (-Inf, u)
    u <- uniroot(f = function(d) pd(d, diff1, sem1, df1, lt = TRUE) - a1,
                 interval = search_int, tol = 1e-06)$root
  } else if (stage == 2) {
    stopifnot(!is.null(diff2), !is.null(sem2), !is.null(df2), !is.null(a2),
              !is.null(weight))
    diff_tmp <- (df1 * diff1 + df2 * diff2) / (df1 + df2)
    search_int <- diff_tmp + c(-6, 6) * max(sem1, sem2)
    # (l, Inf)
    l <- uniroot(f = function(d) comb(pd(d, diff1, sem1, df1, lt = FALSE), 
                                      pd(d, diff2, sem2, df2, lt = FALSE), 
                                      weight) - a2, 
                 interval = search_int, tol = 1e-06)$root
    # (-Inf, u)
    u <- uniroot(f = function(d) comb(pd(d, diff1, sem1, df1, lt = TRUE), 
                                      pd(d, diff2, sem2, df2, lt = TRUE), 
                                      weight) - a2, 
                 interval = search_int, tol = 1e-06)$root
  } else {
    stop("Only stage equal to 1 or 2 implemented.")
  }
  if (l <= u) list(lower = l, upper = u) else NA
}