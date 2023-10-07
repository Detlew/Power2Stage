param.tsd.in <- function(method = c("MCT", "SCT"), alpha = 0.05, 
                         targetpower = 0.8, CV, theta0 = 0.95, theta1, theta2,
                         GMR = theta0, usePE = FALSE, min.n2 = 4, 
                         max.n = list(value = 4000),
                         fCpower = targetpower, fCrit = "CI", fClower, fCupper,
                         fCNmax, ssr.conditional = c("error_power", "error", "no"),
                         nsims = 1E5, setseed = TRUE) {
  
  if (length(alpha) > 1) stop("Length of alpha must be <= 1.")
  method <- match.arg(method)
  ssr.conditional <- match.arg(ssr.conditional)
  if (missing(theta1) && missing(theta2))  theta1 <- 0.8
  if (!missing(theta1) && missing(theta2)) theta2 <- 1/theta1
  if (missing(theta1) && !missing(theta2)) theta1 <- 1/theta2
  
  # Check futility criterion
  fc <- futility_checks_in(fCrit, fClower, fCupper, fCNmax, theta1,
                           n1, min.n2)
  fClower <- fc$fClower
  fCupper <- fc$fCupper
  fCNmax <- fc$fCNmax
  
  # Search for (n1, w) that minimizes the average sample size, given that the
  # overall target power is at least targetpower
  
  # Want to put constraints on n1 and weights. Define them for use in nloptr()
  lw <- switch(method, "SCT" = 1, "MCT" = 2) # length of weights
  # lower bound constraint
  # TO DO: Should we not introduce the limit 24 for CV > 0.3?
  lb <- c(11, rep(0.01, lw))
  n_fixed <- .sampleN2(alpha = alpha, targetpower = targetpower, 
                       ltheta0 = log(theta0), mse = PowerTOST::CV2mse(CV), 
                       method = "nct", bk = 2, dfc = "n-2")
  # upper bound constraint
  ub <- c(max(n_fixed, lb[1] + 2), rep(0.99, lw))
  
  # Now define starting values for optimization algorithm
  x0 <- c((lb[1] + ub[1]) / 2, rep(0.5, lw))
  # We will first perform a global optimization, followed by a local
  # optimization procedure using the starting values from the result of the
  # global algorithm
  alg_g <- "NLOPT_GN_CRS2_LM" # Controlled random search with local mutation
  alg_l <- "NLOPT_LN_NELDERMEAD" # Nelder-Mead Simplex
  
  # Define objective function that needs to be minimized
  # Caveat: For this optimization process we need an allowed maximum sample
  #         size that is finite.
  # We will use max.n = 4000 as upper limit (as used by Maurer et al).
  # Moreover, we will allow the maximum sample size to be a multiple of n1 ("k*n1")
  # or a multiple of the sample size from a fixed study design ("k*nfixed").
  maxn_n1 <- FALSE # is max.n a multiple of n1?
  if (is.null(max.n$type)) {
    stopifnot(is.numeric(max.n$value))
    maxn <- if (is.infinite(max.n$value)) 4000 else max.n$value
    # for later use as return value if power constraint is not met
    max.n.lim <- max.n$value + 1
  } else {
    if (max.n$type == "n1") {
      maxn_n1 <- TRUE
      max.n.lim <- max.n$value * ub[1] + 1
    } else if (max.n$type == "nfixed") {
      maxn <- max.n$value * n_fixed
      max.n.lim <- max.n + 1
    } else {
      stop("Unknown max.n type.")
    }
  }
  
  # Objective function
  f <- function(x, pm = "nct") {
    # function in x = (n1, w)
    n1 <- up2even(x[1])
    if (maxn_n1) maxn <- max.n$value * n1
    w <- x[2:(lw + 1)]
    res <- power.tsd.in(alpha = alpha, weight = w, n1 = n1,
                        max.comb.test = if (lw == 1) FALSE else TRUE,
                        CV = CV, targetpower = targetpower, theta0 = theta0,
                        theta1 = theta1, theta2 = theta2, GMR = GMR,
                        usePE = usePE, min.n2 = min.n2, max.n = maxn,
                        fCpower = fCpower, fCrit = fCrit, fClower = fClower,
                        fCupper = fCupper, fCNmax = fCNmax,
                        ssr.conditional = ssr.conditional, nsims = nsims,
                        setseed = setseed, pmethod = pm)
    if (res$pBE < targetpower - 0.001) max.n.lim else res$nmean
  }
  
  # Perform global optimization
  opts <- list("algorithm" = alg_g, maxeval = 500)
  res <- nloptr::nloptr(x0 = x0, eval_f = f, opts = opts, lb = lb, ub = ub,  
                        pm = "shifted")
  # Set those values as new starting values for local optimization procedure
  opts <- list("algorithm" = alg_l, maxeval = 500)
  res <- nloptr::nloptr(x0 = res$solution, eval_f = f, opts = opts, lb = lb, 
                        ub = ub,  pm = "nct")
  
  if (res$objective == max.n.lim) {
    return(list(n1_opt = NA, w_opt = NA, ass_opt = NA, ratio_n1_nfixed = NA,
                message = NA, status = NA))
  }
  list(n1_opt = up2even(res$solution[1]),
       w_opt = round(res$solution[2:(1 + lw)], 2), 
       ass_opt = res$objective,
       ratio_n1_nfixed = up2even(res$solution[1]) / n_fixed,
       message = res$message, status = res$status)
}