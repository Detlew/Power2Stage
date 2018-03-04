interim.2stage.in <- function(alpha, weight, max.comb.test = TRUE, 
                              targetpower = 0.8, GMR1, n1, CV1, df1 = NULL, 
                              SEM1 = NULL, theta1, theta2, GMR, usePE = FALSE, 
                              min.n2 = 4, max.n = Inf, fCpower = targetpower,
                              fCrit = "CI", fClower, fCupper, fCNmax, 
                              ssr.conditional = c("error_power", "error", "no"),
                              pmethod = c("exact", "nct", "shifted")) {
  
  ### Error handling and default value set-up ----------------------------------
  if (missing(GMR1)) stop("GMR1 must be given.")
  if (missing(CV1))  stop("CV1 must be given!")
  if (CV1 <= 0)      stop("CV1 must be >0!")
  if (missing(n1)) {
    n1 <- NULL
  } else {
    if (n1 <= 0)     stop("Number of subjects in stage 1 must be >0!")
    if (n1 >= max.n) stop("max.n must be greater than n1.")
  }
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
  if (GMR1 <= theta1 || GMR1 >= theta2 || GMR <= theta1 || GMR >= theta2) 
    stop("GMRs must be within acceptance range!")
  
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
  # Test statistics
  t1 <- (lGMR1 - ltheta1) / sem
  t2 <- (lGMR1 - ltheta2) / sem
  # p-values of first stage 1
  p11 <- pt(t1, df = df, lower.tail = FALSE)
  p12 <- pt(t2, df = df, lower.tail = TRUE)
  BE <- (p11 < cl$siglev[1] && p12 < cl$siglev[1])
  
  ## Calculate corresponding exact repeated CI
  # Test inversion for both hypotheses
  f <- function(t)
    pt((lGMR1 - t) / sem, df = df, lower.tail = FALSE) - cl$siglev[1]
  g <- function(t)
    pt((lGMR1 - t) / sem, df = df, lower.tail = TRUE) - cl$siglev[1]
  search_int <- lGMR1 + c(-5, 5) * sem
  ll <- uniroot(f, interval = search_int)$root
  lu <- ll + 2 * (lGMR1 - ll)  # corresponding upper bound
  ru <- uniroot(g, interval = search_int)$root
  rl <- ru - 2 * (ru - lGMR1)  # corresponding lower bound
  # CI for equivalence problem via intersection
  lower <- max(ll, rl)
  upper <- min(lu, ru)
  ci <- if (upper < lower) NA else list(lower = exp(lower), upper = exp(upper))
  
  ## Initialize n2 to be zero (= stop after stage 1)
  n2 <- 0
  
  # Cases with BE == TRUE are clear: early stop due to BE
  # Cases with BE == FALSE are not yet clear:
  # - calculate power for stage 1
  diffm_s1 <- if (usePE) lGMR1 else lGMR
  pwr_s1 <- .calc.power(alpha = cl$siglev[1], ltheta1 = ltheta1, 
                        ltheta2 = ltheta2, diffm = diffm_s1, 
                        sem = sem, df = df, method = pmethod)
  # - if result is FALSE and power for stage 1 is 'sufficiently high', then
  #   the result will be considered a fail (ie leave FALSE)
  #   otherwise (power not high) we leave it open (ie set to NA)
  fut <- vector("integer", 3)
  fut[1] <- (BE == FALSE && pwr_s1 >= fCpower)
  BE[BE == FALSE && pwr_s1 < fCpower] <- NA
  
  # If NA may still consider it as failure due to futility:
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
    fut[2] <- is.na(BE) && outside
    # If identified as futile set to a failure
    BE[is.na(BE) & outside] <- FALSE
    rm(outside)
  }
  
  ### Sample size re-estimation ------------------------------------------------
  if (is.na(BE)) {
    # We need stage 2 (or at least need further information regarding fCNmax)
    
    # Derive components for z-statistics
    Z11 <- qnorm(1 - p11) 
    Z12 <- qnorm(1 - p12)
    
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
    if (fut[3])
      n2 <- 0
  }
  
  ### Define final output ------------------------------------------------------
  res <- list(
    stage = 1, alpha = cl$siglev, cval = cl$cval, weight = weight,
    max.comb.test = max.comb.test, targetpower = targetpower, GMR1 = GMR1, 
    n1 = n1, CV1 = CV1,  df1 = df, SEM1 = sem, theta1 = theta1, theta2 = theta2,
    GMR = GMR, usePE = usePE, min.n2 = min.n2, max.n = max.n, fCpower = fCpower,
    fCrit = fCrit, fCrange = c(fClower, fCupper), fCNmax = fCNmax, 
    ssr.conditional = ssr.conditional, pmethod = pmethod, futility = fut, 
    t11 = t1, t12 = t2, p11 = p11, p12 = p12, 
    CI90 = list(lower = if (nms_match[1]) exp(lower) else NULL,
                upper = if (nms_match[1]) exp(upper) else NULL),
    'Power Stage 1' = pwr_s1,
    n2 = n2, stop_s1 = (n2 == 0), stop_fut = any(fut > 0), 
    RCI = ci, BE = (n2 == 0 && all(fut == 0)), 
    alpha_ssr = if (n2 == 0) NULL else as.numeric(alpha_ssr),
    GMR_ssr = if (n2 == 0) NULL else exp(lGMR_ssr),
    targetpower_ssr = if (n2 == 0) NULL else pwr_ssr
  )
  #class(res) <- c("evaltsd", "list")
  res
}  # end function interim.2stage.in

# alias of function
interim.tsd.in <- interim.2stage.in