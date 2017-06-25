power.2stage.in <- function(alpha = 0.05, n1, CV, targetpower = 0.8, 
                            power.threshold = targetpower, min.n2 = 4, 
                            theta0, theta1, theta2, GMR, #usePE = FALSE,
                            weight = 0.5, max.comb.test = FALSE, 
                            ssr.conditional = TRUE,
                            fCrit = NULL, fClower, fCupper, Nmax, 
                            pmethod = c("nct", "exact", "shifted"), 
                            npct = c(0.05, 0.5, 0.95), nsims, setseed = TRUE, 
                            details = FALSE) {
  # Computes Power or Type I Error rate for the two-stage design scheme
  # based on the inverse normal method. Several design schemes are possible:
  # main scheme is according to Maurer et al
  #
  # Args:
  #   alpha: One-sided significance level
  #   n1: Total sample size for stage 1
  #   CV: Coefficient of variation (use e.g. 0.3 for 30%)
  #   targetpower: Desired target power for end of the trial and power
  #                threshold used for futility check after stage 1
  #   power.threshold: Threshold for power monitoring step to decide on
  #                    futility for cases with not BE after stage 1
  #                    (set to 1 to deactivate this futility rule)
  #   min.n2: Minimum sample size for stage 2
  #   theta0: Assumed ratio of geometric means for simulations
  #   theta1: Lower (bio-)equivalence limits
  #   theta2: Upper (bio-)equivalence limits
  #   GMR: Assumed ratio of geometric means to be used in sample size re-est.
  #   usePE: ...
  #   weight: weight(s) to be used for standard or maximum combination test
  #   max.comb.test: Logical; if TRUE, maximum combination test will be used
  #   ssr.conditional: Logical; if TRUE, the sample size re-estimation step
  #                    uses conditional error rates and conditional power
  #   fCrit: Futility criterion to use: PE, CI or Nmax or a combination thereof
  #   fClower: Lower futility limit for PE or CI of stage 1
  #   fCupper: Upper futility limit for PE or CI of stage 1
  #   Nmax: Maximum allowed total sample size (futility criterion)
  #   pmethod: Power calculation method, to be used in sample size re-est.
  #   npct: Percentiles for the distribution of overall total sample size
  #   nsims: Number of studies to simulate
  #   setseed: Logical; if TRUE, a seed of 1234567 will be set for the sims 
  #   details: logical; if TRUE, print results of time measurements of sims
  #
  # Returns:
  #   Type I Error or Power
  
  # TO DO:
  #  - implement argument usePE
  #  - implement argument max.comb.test
  #  - implement argument details
  #  - implement argument ssr.conditional
  
  ### Error handling and default value set-up ----------------------------------
  alpha <- alpha[1L]
  if (alpha < 0 || alpha > 1) stop("alpha must be within [0, 1]")
  if (missing(n1))  stop("Number of subjects in stage 1 must be given!")
  if (n1 <= 0)      stop("Number of subjects in stage 1 must be >0!")
  if (missing(CV))  stop("CV must be given!")
  if (CV <= 0)      stop("CV must be >0!")
  if (targetpower < 0 || targetpower > 1) 
    stop("targetpower must be within [0, 1]")
  if (power.threshold < 0 || power.threshold > 1) 
    stop("power_threshold must be within [0, 1]")
  if (min.n2 < 4) 
    stop("min.n2 has to be at least 4.")
  if (min.n2 %% 2 != 0) {  # make it even
    min.n2 <- min.n2 + min.n2 %% 2
    message("min.n2 rounded up to next even integer", min.n2)
  }
  if (missing(GMR)) GMR <- 0.95
  if (missing(theta0)) theta0 <- GMR
  if (missing(theta1) && missing(theta2))  theta1 <- 0.8
  if (!missing(theta1) && missing(theta2)) theta2 <- 1/theta1
  if (missing(theta1) && !missing(theta2)) theta1 <- 1/theta2
  if (GMR <= theta1 || GMR >= theta2) 
    stop("GMR must be within acceptance range!")
  
  # TO DO:
  #  - add check for argument weight
  #  - extend it for weight being of length 2 (for max.comb.test == TRUE)
  
  # Check futility criterion
  if (!is.null(fCrit)) {
    fCrit <- tolower(fCrit)
    fcrit_nms <- c("pe", "ci", "nmax")  # correct possibilities
    nms_match <- fcrit_nms %in% fCrit  # check which fCrits are given
    if (sum(nms_match) == 0)
      stop("fCrit not correctly specified.")
    if (nms_match[1]) {  # PE
      if (missing(fClower) && missing(fCupper))  fClower <- theta1
      if (missing(fClower) && !missing(fCupper)) fClower <- 1/fCupper
      if (!missing(fClower) && missing(fCupper)) fCupper <- 1/fClower
    }
    if (nms_match[2]) {  # CI
      if (nms_match[1]) {
        message("Both PE and CI specified for futility, PE will be ignored.")
        nms_match[1] <- FALSE
      }
      if (missing(fClower) && missing(fCupper))  fClower <- 0.95
      if (missing(fClower) && !missing(fCupper)) fClower <- 1/fCupper
      if (!missing(fClower) && missing(fCupper)) fCupper <- 1/fClower
    }
    if (nms_match[3]) {  # Nmax
      if (missing(Nmax)) Nmax <- 4*n1
      if (!missing(Nmax) && (Nmax < n1 + min.n2))
        stop("Nmax must be greater than n1 + min.n2.")
      if (sum(nms_match[1:2]) == 0) {
        fClower <- 0
        fCupper <- Inf
      }
    } else {
      Nmax <- Inf
    }
  } else {
    fClower <- 0
    fCupper <- Inf
    Nmax <- Inf
  }
  
  # Check if power calculation method is nct, exact or shifted
  pmethod <- match.arg(pmethod)
  # Set default value for nsims if missing
  if (missing(nsims))
    nsims <- if (theta0 <= theta1 || theta0 >= theta2) 1E6 else 1E5
  if (setseed) set.seed(1234567)
  ltheta0 <- log(theta0)
  ltheta1 <- log(theta1)
  ltheta2 <- log(theta2)
  lGMR    <- log(GMR)
  mse     <- CV2mse(CV)
  bk      <- 2  # 2x2x2 crossover design const
  ntot    <- rep.int(n1, nsims)  # initialize overall sample size
  
  ### Calculate adjusted significance level ------------------------------------
  # TO DO:
  #  - to be adapted when weight is more than length 1 
  #bnds <- ldbounds::bounds(t = c(weight, 1), iuse = 2, alpha = alpha)
  #alpha_adj <- 1 - pnorm(bnds$upper.bounds)
  cvals <- critical.value.2stage(alpha, weight)
  
  ### Stage 1 ------------------------------------------------------------------
  stage <- rep.int(1, nsims)
  
  ## Simulation of stage 1 data
  se.fac <- sqrt(bk / n1)
  df <- n1 - 2
  
  # Simulate point estimates via normal distribution
  pes <- rnorm(n = nsims, mean = ltheta0, sd = se.fac * sqrt(mse))
  # and MSE via chi-squared distribution
  mses <- mse * rchisq(n = nsims, df = df) / df
  
  # The 2 one-sided test statistics
  t1 <- (pes - ltheta1) / sqrt(mses) / se.fac
  t2 <- (pes - ltheta2) / sqrt(mses) / se.fac
  
  # p-values of first stage
  p11 <- pt(t1, df = df, lower.tail = FALSE)
  p12 <- pt(t2, df = df, lower.tail = TRUE)
  
  rm(t1, t2)
  
  ## Evaluation of stage 1
  BE <- (p11 < cvals$siglev & p12 < cvals$siglev)
  
  # Cases with BE == TRUE are clear: early stop due to BE
  # Cases with BE == FALSE are not yet clear:
  # - calculate power for stage 1
  pwr_s1 <- .calc.power(alpha = cvals$siglev, ltheta1 = ltheta1, 
                        ltheta2 = ltheta2, diffm = lGMR, 
                        sem = se.fac * sqrt(mses), df = df, method = pmethod)
  # - if result is FALSE and power for stage 1 is 'sufficiently high', then
  #   the result will be considered a fail (ie leave FALSE)
  #   otherwise we leave it open (ie set to NA)
  BE[BE == FALSE & pwr_s1 < power.threshold] <- NA
  
  # From those NAs may still consider some of them as failure due to futility:
  # Futility check - here only regarding PE or CI (Nmax comes later)
  if (!is.null(fCrit) && sum(nms_match[1:2]) > 0) {
    lfClower <- log(fClower)
    lfCupper <- log(fCupper)
    if (nms_match[1]) {
      outside <- ((pes - lfClower) < 1.25e-5 | (lfCupper - pes) < 1.25e-5)
    }
    if (nms_match[2]) {
      tval <- qt(1 - 0.05, df)  # use 90% CI for this check
      hw <- tval * se.fac * sqrt(mses)
      lower <- pes - hw
      upper <- pes + hw
      outside <- (lower > lfCupper) | (upper < lfClower)
      rm(tval, hw, lower, upper)
    }
    # Set the ones identified as futile to a failure
    BE[is.na(BE) & outside] <- FALSE
    rm(outside)
  }
  
  ### Sample size re-estimation ------------------------------------------------
  # Only consider scenarios where stage 2 is necessary
  pes_tmp  <- pes[is.na(BE)]
  mses_tmp <- mses[is.na(BE)]
  pwr_s1   <- pwr_s1[is.na(BE)]
  p11      <- p11[is.na(BE)]
  p12      <- p12[is.na(BE)]
  
  # Do it only if stage 2 is required for at least one scenario
  if (length(pes_tmp) > 0) {
    # Set stage variable to 2 for the relevant scenarios
    # (NB: some may be changed to 1 again after all -> Nmax criterion)
    stage[is.na(BE)] <- 2
    
    # Define temporary variables for BE and stage for the second stage
    BE2 <- rep.int(NA, length(pes_tmp))
    s2 <- rep.int(2, length(pes_tmp))
    
    # Derive components for z-statistics
    Z11 <- qnorm(1 - p11) 
    Z12 <- qnorm(1 - p12)
    
    if (ssr.conditional) {
      # Calculate conditional power (may be negative, if so, set to zero)
      pwr_ssr <- pmax(((1 - pwr_s1) - (1 - targetpower)) / (1 - pwr_s1), 0)
    
      # Derive conditional error rates
      # TO DO: 
      #  - adapt alpha_c in case of max.comb.test == TRUE
      alpha_ssr <- 1 - pnorm((cvals$cval - sqrt(weight)*cbind(Z11, Z12)) / 
                             sqrt(1 - weight))
    
      # Set sign of lGMR to the sign of estimated point estimate
      # (Maurer et al call this 'adaptive planning step')
      sgn_pes_tmp <- ifelse(pes_tmp >= 0, 1, -1)
      lGMR_ssr <- lGMR * sgn_pes_tmp * if (lGMR >= 0) 1 else -1
    } else {
      alpha_ssr <- cvals$siglev  # alpha_adj[2]
      pwr_ssr <- targetpower
      lGMR_ssr <- lGMR  # here we do not need to switch the sign (?!)
    }
    
    # Sample size for stage 2
    n2 <- .sampleN3(alpha = alpha_ssr, targetpower = pwr_ssr, 
                    ltheta0 = lGMR_ssr, mse = mses_tmp, 
                    ltheta1 = ltheta1, ltheta2 = ltheta2, method = pmethod)
    if (!ssr.conditional)
      n2 <- n2 - n1
    n2 <- pmax(n2, min.n2)
    
    # Futility check regarding maximum overall sample size
    # - if n2 is infinite, then we consider this as fail too
    BE2[is.infinite(n2)] <- FALSE
    # - Otherwise check if n1+n2 is greater than Nmax (if so, set to fail)
    if (!is.null(fCrit) && nms_match[3])
      BE2[n1 + n2 > Nmax] <- FALSE
    s2[BE2 == FALSE] <- 1  # such a case is considered to be gone up to stage 1
    
    # Carry over results from BE2 and s2 to BE, stage and ntot
    stage[is.na(BE)] <- s2
    ntot[is.na(BE)] <- n1 + n2
    BE[is.na(BE)] <- BE2
    
    # Reduce relevant characteristics to the ones where stage 2 will be done
    p11 <- p11[is.na(BE2)]
    p12 <- p12[is.na(BE2)]
    Z11 <- Z11[is.na(BE2)]
    Z12 <- Z12[is.na(BE2)]
    n2  <- n2[is.na(BE2)]
    
    rm(BE2, s2, alpha_ssr, pes_tmp, mses_tmp, pwr_ssr, lGMR_ssr, pwr_s1)
    
    ### Stage 2 ----------------------------------------------------------------
    ## Simulation of stage 2 data
    se.fac <- sqrt(bk / n2)
    df <- n2 - 2
    nsims2 <- length(p11)
    
    # Simulate point estimates via normal distribution
    pes_s2 <- rnorm(n = nsims2, mean = ltheta0, sd = se.fac * sqrt(mse))
    # and MSE via chi-squared distribution
    mses_s2 <- mse * rchisq(n = nsims2, df = df) / df
    
    # The 2 one-sided test statistics
    t1 <- (pes_s2 - ltheta1) / sqrt(mses_s2) / se.fac
    t2 <- (pes_s2 - ltheta2) / sqrt(mses_s2) / se.fac
    
    # p-values of second stage
    p21 <- pt(t1, df = df, lower.tail = FALSE)
    p22 <- pt(t2, df = df, lower.tail = TRUE)
    
    rm(pes_s2, mses_s2, t1, t2)
   
    ## Combined evaluation via inverse normal approach
    # Derive further components for z-statistics 
    Z21 <- qnorm(1 - p21)
    Z22 <- qnorm(1 - p22)
    # For H01, test statistic at the end of second stage:
    Z01 <- sqrt(weight) * Z11 + sqrt(1 - weight) * Z21
    # For H02, test statistic at the end of second stage:
    Z02 <- sqrt(weight) * Z12 + sqrt(1 - weight) * Z22
    
    # Fill the remaining NAs with either TRUE or FALSE
    # TO DO: max.comb.test boundary
    BE[is.na(BE)] <- (Z01 > cvals$cval & Z02 > cvals$cval)
    
    rm(Z11, Z21, Z12, Z22, Z01, Z02)
  }  # end if length(pes_tmp) > 0
  
  ### Define final output ------------------------------------------------------
  res <- list(
    # Information
    design = "2x2 crossover", alpha = cvals$siglev, n1 = n1, CV = CV, GMR = GMR, 
    targetpower = targetpower, min.n2 = min.n2, theta0 = theta0, 
    theta1 = theta1, theta2 = theta2, weight = weight, 
    max.comb.test = max.comb.test, fCrit = fCrit, fCrange = c(fClower, fCupper),
    Nmax = Nmax, pmethod = pmethod, nsims = nsims,
    # Results
    pBE = sum(BE) / nsims,
    pBE_s1 = sum(BE[ntot == n1]) / nsims,
    pct_s2 = 100 * sum(ntot > n1) / nsims,
    pct_stop_s1 = 100 * sum(ntot == n1) / nsims,
    nmean = mean(ntot),
    nrange = range(ntot),
    nperc = quantile(ntot, probs = npct)
  )
  # TO DO: implement class attributes properly:
  #  - include pct_stop_s1 in order to get stopping rates as in Maurer et al
  #  - output of futility criterion to be adapted
  #  - output of n(s1, s2) to be adapted
  #class(res) <- c("pwrtsd", "list")
  res
}  # end function power.2stage.in