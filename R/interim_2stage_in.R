interim.2stage.in <- function(alpha, weight, max.comb.test = TRUE, 
                              targetpower = 0.8, GMR1, n1, CV1, df1 = NULL, 
                              SEM1 = NULL, theta1, theta2, GMR, usePE = FALSE, 
                              min.n2 = 4, max.n = Inf, fCpower = targetpower,
                              fCrit = "CI", fClower, fCupper, fCNmax, 
                              ssr.conditional = c("error_power", "error", "no"),
                              pmethod = c("exact", "nct", "shifted")) {
  
  ### Error handling and default value set-up ----------------------------------
  if (missing(GMR1)) stop("GMR1 must be given.")
  if (missing(CV1))  stop("CV1 must be given.")
  if (CV1 <= 0)      stop("CV1 must be >0.")
  if (missing(n1))   stop("n1 must be given.")
  if (n1 <= 0)       stop("Number of subjects in stage 1 must be >0.")
  if (n1 >= max.n)   stop("max.n must be greater than n1.")
  if (missing(alpha))
    alpha <- 0.05
  if (length(alpha) > 2) 
    stop("Length of alpha must be <= 2.")
  if (!missing(weight) && length(weight) > 2)
    stop("Length of weight must be <= 2.")
  if (length(alpha) == 1 && missing(weight))
    weight <- if (max.comb.test) c(0.5, 0.25) else 0.5
  if (length(alpha) == 2) {
    if (missing(weight)) {
      stop("weight must be specified.")
    } else {
      message(paste0("Note: Adjusted alphas are specified, it is assumed that",
                     " the specified weight(s) are in line with the alpha",
                     " values."))
    }
  }
  lw <- length(weight)
  if (max.comb.test) {
    if (lw != 2)
      stop("Two weights, w and w*, are required for maximum combination test.")
  } else {
    if (lw != 1)
      stop("One weight, w, is required for standard combination test.")
  }
  if (any(weight <= 0) || any(weight >= 1))
    stop("Weight(s) not properly specified, must be > 0 and < 1.")
  if (any(alpha <= 0) || any(alpha >= 1))
    stop("Alpha(s) not properly specified, must be > 0 and < 1.")
  if (targetpower <= 0 || targetpower >= 1) 
    stop("targetpower must be within (0, 1)")
  if (fCpower < 0 || fCpower > 1) 
    stop("fCpower must be within [0, 1]")
  if (min.n2 < 4) 
    stop("min.n2 has to be at least 4.")
  if (min.n2 %% 2 != 0) {  # make it even
    min.n2 <- 2 * floor(min.n2 / 2) + 2
    message("min.n2 rounded up to next even integer", min.n2)
  }
  if (is.finite(max.n) && (max.n %% 2 != 0)) {
    max.n <- 2 * floor(max.n / 2) + 2
    message("max.n rounded up to next even integer", max.n)
  }
  if (missing(GMR)) GMR <- 0.95
  if (missing(theta1) && missing(theta2))  theta1 <- 0.8
  if (!missing(theta1) && missing(theta2)) theta2 <- 1/theta1
  if (missing(theta1) && !missing(theta2)) theta1 <- 1/theta2
  
  # Check futility criterion
  stopifnot(is.character(fCrit))
  fCrit <- tolower(fCrit)
  fcrit_nms <- c("ci", "pe", "nmax", "no")  # correct possibilities
  nms_match <- fcrit_nms %in% fCrit  # check which fCrits are given
  if (sum(nms_match) == 0)
    stop("fCrit not correctly specified.")
  if (nms_match[4]) { # No futility criterion
    if (sum(nms_match[1:3]) > 0) {
      message("No futility will be applield.")
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
  ssr.conditional <- match.arg(ssr.conditional)
  pmethod <- match.arg(pmethod)
  lGMR1   <- log(GMR1)
  lGMR    <- log(GMR)
  ltheta1 <- log(theta1)
  ltheta2 <- log(theta2)
  mse     <- CV2mse(CV1)
  des <- get_n_df_sem(n = n1, df = df1, mse = mse, sem = SEM1)
  n1 <- des$n
  df <- des$df
  sem <- des$sem
  
  ### Calculate adjusted critical levels ---------------------------------------
  cl <- if (length(alpha) == 1) critical.value.2stage(alpha, weight) else
    list(cval = qnorm(1 - alpha), siglev = alpha)
  
  ### Evaluation of Stage 1 ----------------------------------------------------

  ## Bioequivalence?
  t1 <- (lGMR1 - ltheta1) / sem
  t2 <- (lGMR1 - ltheta2) / sem
  # p-values of first stage 1
  p11 <- pt(t1, df = df, lower.tail = FALSE)
  p12 <- pt(t2, df = df, lower.tail = TRUE)
  # Derive components for z-statistics (for later use)
  Z11 <- qnorm(1 - p11) 
  Z12 <- qnorm(1 - p12)
  # BE yes/no
  BE <- (p11 <= cl$siglev[1] && p12 <= cl$siglev[1])
  
  ## Calculate corresponding exact repeated CI
  rci <- repeated_ci(diff1 = lGMR1, sem1 = sem, df1 = df, a1 = cl$siglev[1],
                     stage = 1)
  
  # Calculate power for stage 1
  diffm_s1 <- lGMR # if (usePE) lGMR1 else lGMR
  pwr_s1 <- .calc.power(alpha = cl$siglev[1], ltheta1 = ltheta1, 
                        ltheta2 = ltheta2, diffm = diffm_s1, 
                        sem = sem, df = df, method = pmethod)
  fut <- vector("integer", 3)
  fut[1] <- (BE == FALSE && pwr_s1 >= fCpower)
  
  # Futility check - here only regarding PE or CI (fCNmax comes later)
  if (fCrit != "no" && sum(nms_match[1:2]) > 0) {
    lfClower <- log(fClower)
    lfCupper <- log(fCupper)
    if (nms_match[1]) {
      tval <- qt(1 - 0.05, df)  # use 90% CI
      hw <- tval * sem
      lower <- lGMR1 - hw
      upper <- lGMR1 + hw
      outside <- (lower > lfCupper) || (upper < lfClower)
    }
    if (nms_match[2]) {
      outside <- ((lGMR1 - lfClower) < 1.25e-5 | (lfCupper - lGMR1) < 1.25e-5)
    }
    fut[2] <- outside
    rm(outside)
  }
  
  n2 <- 0
  if (!BE) {
    ### Sample size re-estimation ------------------------------------------------
    if (ssr.conditional == "no") {
      alpha_ssr <- cl$siglev[2]
      pwr_ssr <- targetpower
      lGMR_ssr <- if (usePE) lGMR1 else lGMR
    } else {
      # Derive conditional error rates
      alpha_ssr <- 1 - pnorm(pmin(
        (cl$cval[2] - sqrt(weight[1])*cbind(Z11, Z12)) / sqrt(1 - weight[1]),
        (cl$cval[2] - sqrt(weight[lw])*cbind(Z11, Z12)) / sqrt(1 - weight[lw])
      ))
    
      # Define target power for ssr
      pwr_ssr <- targetpower
      if ((ssr.conditional == "error_power") && (fCpower <= targetpower)) {
        # Use conditional power
        pwr_ssr <- 1 - (1 - targetpower) / (1 - pwr_s1)
      }
    
      if (usePE) {
        lGMR_ssr <- lGMR1
        if (GMR1 <= theta1 || GMR1 >= theta2) {
          message(paste0("SSR using observed GMR being outside of ", 
                         "theta1 ... theta2 not possible, use planned GMR."))
          lGMR_ssr <- lGMR
        }
      } else {
        # Set sign of lGMR to the sign of estimated point estimate
        # (Maurer et al call this 'adaptive planning step')
        sgn_pe <- if (lGMR1 >= 0) 1 else -1
        lGMR_ssr <- abs(lGMR) * sgn_pe
      }
    }
  
    # Sample size for stage 2
    n2 <- .sampleN3(alpha = alpha_ssr, targetpower = pwr_ssr, 
                    ltheta0 = lGMR_ssr, mse = mse, 
                    ltheta1 = ltheta1, ltheta2 = ltheta2, method = pmethod)
    if (ssr.conditional == "no")
      n2 <- n2 - n1
    n2 <- max(min(n2, max.n - n1), min.n2)
    
    # Futility check regarding maximum overall sample size
    fut[3] <- (n1 + n2 > fCNmax) || is.infinite(n2)
  }
  
  ### Define final output ------------------------------------------------------
  ci90 <- if (nms_match[1]) exp(c(lower, upper)) else NULL
  if (!is.null(ci90)) names(ci90) <- c("lower CL", "upper CL")
  res <- list(
    stage = 1L, alpha = cl$siglev, cval = cl$cval, weight = weight,
    max.comb.test = max.comb.test, targetpower = targetpower, GMR1 = GMR1, 
    n1 = as.integer(n1), CV1 = CV1, df1 = df, SEM1 = sem, theta1 = theta1, 
    theta2 = theta2, GMR = GMR, usePE = usePE, min.n2 = as.integer(min.n2),
    max.n = if (is.infinite(max.n)) Inf else as.integer(max.n), 
    fCpower = fCpower, fCrit = fCrit, fCrange = c(fClower, fCupper), 
    fCNmax = if (is.infinite(fCNmax)) Inf else as.integer(fCNmax),
    ssr.conditional = ssr.conditional, pmethod = pmethod,
    #t11 = t1, t12 = t2,
    p11 = p11, p12 = p12, z1 = Z11, z2 = Z12,
    futility = fut,
    CI90 = ci90,
    'Power Stage 1' = pwr_s1, n2 = as.integer(n2), 
    stop_s1 = (BE == TRUE) || any(fut > 0),
    stop_fut = any(fut > 0), stop_BE = (BE == TRUE), RCI = exp(rci), 
    alpha_ssr = if (!BE) as.numeric(alpha_ssr) else NULL, 
    GMR_ssr = if (!BE) exp(lGMR_ssr) else NULL,
    targetpower_ssr = if (!BE) pwr_ssr else NULL
  )
  class(res) <- c("evaltsd", "list")
  res
}  # end function interim.2stage.in

# alias of function
interim.tsd.in <- interim.2stage.in