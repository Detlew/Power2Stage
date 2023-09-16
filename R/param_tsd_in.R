up2even <- function(n) {
  if (n %% 2 != 0) return(2 * ceiling(n/2))
  n
}

param.tsd.in <- function(method = c("SCT", "MCT"), alpha = 0.05, 
                         targetpower = 0.8, CV, theta0 = 0.95, theta1, theta2,
                         GMR = theta0, usePE = FALSE, min.n2 = 4, max.n = 4000,
                         fCpower = targetpower, fCrit = "CI", fClower, fCupper,
                         fCNmax, ssr.conditional = c("error_power", "error", "no"),
                         nsims = 1E5, setseed = TRUE) {
  
  if (length(alpha) > 1) stop("Length of alpha must be <= 1.")
  method <- match.arg(method)
  ssr.conditional <- match.arg(ssr.conditional)
  if (missing(theta1) && missing(theta2))  theta1 <- 0.8
  if (!missing(theta1) && missing(theta2)) theta2 <- 1/theta1
  if (missing(theta1) && !missing(theta2)) theta1 <- 1/theta2
  
  # TO DO: Without the following checks on futility criterion the function
  #        call later won't work -> check if this can be avoided somehow as it
  #        is lengthy and already part of power.tsd.in()
  stopifnot(is.character(fCrit))
  fCrit <- tolower(fCrit)
  fcrit_nms <- c("ci", "pe", "nmax", "no")  # correct possibilities
  nms_match <- fcrit_nms %in% fCrit  # check which fCrits are given
  if (sum(nms_match) == 0)
    stop("fCrit not correctly specified.")
  if (nms_match[4]) { # No futility criterion
    if (sum(nms_match[1:3]) > 0) {
      message("No futility will be applied.")
    }
    fClower <- 0
    fCupper <- Inf
    fCNmax <- Inf
  } else {
    if (nms_match[1]) {  # CI
      if (nms_match[2]) {
        message("Both PE and CI specified for futility, PE will be ignored.")
        nms_match[1] <- FALSE
      }
      if (missing(fClower) && missing(fCupper))  fClower <- 0.95
      if (missing(fClower) && !missing(fCupper)) fClower <- 1/fCupper
      if (!missing(fClower) && missing(fCupper)) fCupper <- 1/fClower
    }
    if (nms_match[2]) {  # PE
      if (missing(fClower) && missing(fCupper))  fClower <- theta1
      if (missing(fClower) && !missing(fCupper)) fClower <- 1/fCupper
      if (!missing(fClower) && missing(fCupper)) fCupper <- 1/fClower
    }
    if (nms_match[3]) {  # Nmax
      if (missing(fCNmax)) fCNmax <- 4*n1
      if (!missing(fCNmax) && (fCNmax < n1 + min.n2))
        stop("fCNmax must be greater than n1 + min.n2.")
      if (sum(nms_match[1:2]) == 0) {
        fClower <- 0
        fCupper <- Inf
      }
    } else {
      fCNmax <- Inf
    }
  }
  
  # Search for (n1, w) that minimize average sample size, given that the
  # overall target power is at least targetpower
  
  # Want constraints for n1 and weights, lb and ub respectively
  # define them for use in nloptr()
  lw <- switch(method, "SCT" = 1, "MCT" = 2) # length of weights
  # TO DO: Should we not introduce the limit 24 for CV > 0.3?
  lb <- c(-1 + if (CV <= 0.3) 12 else 24, rep(0.01, lw))
  n_fixed <- .sampleN2(alpha = alpha, targetpower = targetpower, 
                       ltheta0 = log(theta0), mse = PowerTOST::CV2mse(CV), 
                       method = "nct", bk = 2, dfc = "n-2")
  ub <- c(max(n_fixed, lb[1] + 2), rep(0.99, lw))
  # Now define starting values
  x0 <- c((lb[1] + ub[1]) / 2, rep(0.5, lw))
  # We will first perform a global optimization, followed by a local
  # optimization procedure usingn the starting values from the result of the
  # global algorithm
  alg_g <- "NLOPT_GN_CRS2_LM" # Controlled random search with local mutation
  alg_l <- "NLOPT_LN_NELDERMEAD" # Nelder-Mead Simplex
  
  # Define objective function that needs to be minimized
  # Caveat: For this optimization process we need an allowed maximum sample
  #         size that is finite.
  
  # We will use max.n = 4000 as upper limit (as used by Maurer et al).
  # Moreover, we will allow the maximum sample size to be a multiple of n1 ("k*n1")
  # or a multiple of the sample size from a fixed study design ("k*nfixed").
  # TO DO: Can this be made more elegant (than is.numeric and is.character)?
  if (is.numeric(max.n)) {
    if (is.infinite(max.n)) max.n <- 4000
    # for later use as return value if power constraint is not met
    max.n.lim <- max.n + 1
  }
  maxn_n1 <- FALSE  # should max n be a multiple of n1?
  if (is.character(max.n)) {
    # Only "k*n1" or "k*nfixed" allowed
    k <- as.numeric(substr(max.n, 1, 1))
    maxn_txt <- substr(max.n, 3, nchar(max.n))
    if (maxn_txt == "nfixed") {
      max.n <- k * n_fixed
      max.n.lim <- max.n + 1
    } else if (maxn_txt == "n1") {
      maxn_n1 <- TRUE
      max.n.lim <- k * ub[1] + 1
    } else {
      stop("Unknown max.n phrase.")
    }
  }
  
  # Objective function
  f <- function(x, pm = "nct") {
    # function in x = (n1, w)
    n1 <- up2even(x[1])
    if (maxn_n1) max.n <- k * n1
    w <- x[2:(lw + 1)]
    res <- power.tsd.in(alpha = alpha, weight = w, n1 = n1,
                        max.comb.test = if (lw == 1) FALSE else TRUE,
                        CV = CV, targetpower = targetpower, theta0 = theta0,
                        theta1 = theta1, theta2 = theta2, GMR = GMR,
                        usePE = usePE, min.n2 = min.n2, max.n = max.n,
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