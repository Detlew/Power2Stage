################################################################################
# Author: Benjamin Lang
# Code based on power.2stage() and power.2stage.fC() by Detlew Labes
################################################################################
power.2stage.in <- function(alpha, weight, max.comb.test = TRUE, n1, CV, 
                            targetpower = 0.8, theta0, theta1, theta2, 
                            GMR, usePE = FALSE, min.n2 = 4, max.n = Inf, 
                            fCpower = targetpower, 
                            fCrit = "CI", fClower, fCupper, fCNmax, 
                            ssr.conditional = c("error_power", "error", "no"),
                            pmethod = c("nct", "exact", "shifted"), 
                            npct = c(0.05, 0.5, 0.95), nsims, setseed = TRUE, 
                            details = FALSE) {
  # Computes Power or Type I Error rate for the two-stage design scheme
  # based on the inverse normal method. Several design schemes are possible:
  # main scheme is according to Maurer et al
  #
  # Args:
  #   alpha: Overall one-sided significance level, if of length 1
  #          Adjusted one-sided alpha levels for stage 1 and 2, if of length 2
  #   weight: Required if alpha is of length 1: weight of first stage
  #           Of length 1 in case of max.comb.test = FALSE, length 2 otherwise
  #   max.comb.test: Logical; if TRUE, maximum combination test will be used
  #   n1: Sample size for stage 1
  #   CV: Coefficient of variation (use e.g. 0.3 for 30%)
  #   targetpower: Desired target power for end of the trial 
  #   theta0: Assumed ratio of geometric means for simulations
  #   theta1: Lower (bio-)equivalence limits
  #   theta2: Upper (bio-)equivalence limits
  #   GMR: Assumed ratio of geometric means to be used in sample size re-est.
  #   usePE: Logical; if TRUE, ssr is done with observed point estimate from
  #          stage 1
  #   min.n2: Minimum sample size for stage 2
  #   max.n: Maximum overall sample size (stage 1 + stage 2) 
  #   fCpower: Threshold for power monitoring step to decide on futility 
  #            for cases 'not BE' after stage 1
  #   fCrit: Futility criterion to use: CI, PE, Nmax or No
  #   fClower: Lower futility limit for PE or CI of stage 1
  #   fCupper: Upper futility limit for PE or CI of stage 1
  #   fCNmax: If re-estimated sample size n2 is such that n1+n2 > fCNmax,
  #           the study will not continue to stage 2 and will be considered
  #           an overall failure (no BE)
  #   ssr.conditional: Specifies the use of conditional error rates and
  #                    conditional power, only conditional error rates or
  #                    no use at all of conditional values in the sample size 
  #                    re-estimation step after stage 1
  #   pmethod: Power calculation method to be used in sample size re-est.
  #   npct: Percentiles for the distribution of overall (total) sample size
  #   nsims: Number of studies to simulate
  #   setseed: Logical; if TRUE, a seed of 1234567 will be set for the sims 
  #   details: logical; if TRUE, print results of time measurements of sims
  #
  # Returns:
  #   Type I Error or Power
  
  ### Error handling and default value set-up ----------------------------------
  if (missing(alpha)) 
    alpha <- 0.05
  if (length(alpha) > 2) 
    stop("Length of alpha must be <= 2.")
  if (!missing(weight) && length(weight) > 2)
    stop("Length of weight must be <= 2.")
  if (length(alpha) == 1 && missing(weight))
    weight <- if (max.comb.test) c(0.5, 0.25) else 0.5
  if (length(alpha) == 2) {
    if (missing(weight))
      stop("weight must be specified.")
    message(paste0("Note: It is assumed that the specified alpha values are in", 
                   " line with the specified max.comb.test argument."))
  }
  if (length(alpha) == 2 && !missing(weight))
    message(paste0("Note: Adjusted alphas are specified, it is assumed that the",
                   " specified weight(s) are in line with the alpha values."))
  lw <- length(weight)
  if (max.comb.test) {
    if (length(alpha) != 2 && lw != 2)
      stop("Two weights, w and w*, are required for maximum combination test.")
  } else {
    if (length(alpha) == 1 && lw != 1)
      stop("One weight, w, is required for standard combination test.")
  }
  if (any(weight <= 0) || any(weight >= 1))
    stop("Weight(s) not properly specified, must be > 0 and < 1.")
  if (any(alpha <= 0) || any(alpha >= 1))
    stop("Alpha(s) not properly specified, must be > 0 and < 1.")
  if (missing(n1))  stop("Number of subjects in stage 1 must be given!")
  if (n1 <= 0)      stop("Number of subjects in stage 1 must be >0!")
  if (n1 %% 2 != 0) {  # make it even
    n1 <- 2 * floor(n1 / 2) + 2
    message("n1 rounded up to next even integer ", n1)
  }
  if (missing(CV))  stop("CV must be given!")
  if (CV <= 0)      stop("CV must be >0!")
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
  if (n1 >= max.n) stop("max.n must be greater than n1.")
  if (is.finite(max.n) && (max.n %% 2 != 0)) {
    max.n <- 2 * floor(max.n / 2) + 2
    message("max.n rounded up to next even integer", max.n)
  }
  if (missing(GMR)) GMR <- 0.95
  if (missing(theta0)) theta0 <- GMR
  if (missing(theta1) && missing(theta2))  theta1 <- 0.8
  if (!missing(theta1) && missing(theta2)) theta2 <- 1/theta1
  if (missing(theta1) && !missing(theta2)) theta1 <- 1/theta2
  if (GMR <= theta1 || GMR >= theta2) 
    stop("GMR must be within acceptance range!")
  
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
  
  # Check usage of conditional error rates / power
  ssr.conditional <- match.arg(ssr.conditional)
  
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
  
  ### Calculate adjusted critical levels ---------------------------------------
  # length(alpha)=1 => if length(weight)=2, max.comb.test=TRUE is assumed here
  cl <- if (length(alpha) == 1) critical.value.2stage(alpha, weight) else
    # if two alphas, ignore possibly given weights and only use those alphas
    list(cval = qnorm(1 - alpha), siglev = alpha)
  
  ### Stage 1 ------------------------------------------------------------------
  if (details) {
    cat(nsims,"sims. Stage 1")
    # Start timer
    ptm  <- proc.time()
  }
  
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
  BE <- (p11 < cl$siglev[1] & p12 < cl$siglev[1])
  
  # Cases with BE == TRUE are clear: early stop due to BE
  # Cases with BE == FALSE are not yet clear:
  # - calculate power for stage 1
  diffm_s1 <- if (usePE) pes else lGMR
  pwr_s1 <- .calc.power(alpha = cl$siglev[1], ltheta1 = ltheta1, 
                        ltheta2 = ltheta2, diffm = diffm_s1, 
                        sem = se.fac * sqrt(mses), df = df, method = pmethod)
  # - if result is FALSE and power for stage 1 is 'sufficiently high', then
  #   the result will be considered a fail (ie leave FALSE)
  #   otherwise (power not high) we leave it open (ie set to NA)
  fut <- sum(BE == FALSE & pwr_s1 >= fCpower)
  BE[BE == FALSE & pwr_s1 < fCpower] <- NA
  
  # From those NAs may still consider some of them as failure due to futility:
  # Futility check - here only regarding PE or CI (fCNmax comes later)
  if (fCrit != "no" && sum(nms_match[1:2]) > 0) {
    lfClower <- log(fClower)
    lfCupper <- log(fCupper)
    if (nms_match[1]) {
      tval <- qt(1 - 0.05, df)  # use 90% CI
      #tval <- qt(1 - cl$siglev[1], df)  # use adjusted CI for this check
      hw <- tval * se.fac * sqrt(mses)
      lower <- pes - hw
      upper <- pes + hw
      outside <- (lower > lfCupper) | (upper < lfClower)
      rm(tval, hw, lower, upper)
    }
    if (nms_match[2]) {
      outside <- ((pes - lfClower) < 1.25e-5 | (lfCupper - pes) < 1.25e-5)
    }
    fut <- fut + sum(is.na(BE) & outside)
    # Set the ones identified as futile to a failure
    BE[is.na(BE) & outside] <- FALSE
    rm(outside)
  }
  
  # Time for stage 1
  if (details) {
    cat(" - Time consumed (secs):\n")
    print(round((proc.time() - ptm), 1))
  }
  
  ### Sample size re-estimation ------------------------------------------------
  # Only consider scenarios where stage 2 is necessary
  pes_tmp  <- pes[is.na(BE)]
  mses_tmp <- mses[is.na(BE)]
  pwr_s1   <- pwr_s1[is.na(BE)]
  p11      <- p11[is.na(BE)]
  p12      <- p12[is.na(BE)]
  rm(pes, mses)
  
  # Do it only if stage 2 is required for at least one scenario
  if (length(pes_tmp) > 0) {
    if (details) {
      cat("Keep calm. Sample sizes for stage 2 (", length(pes_tmp),
          " studies)\n", sep = "")
      cat("will be estimated. May need some time.\n")
    }
    
    # Set stage variable to 2 for the relevant scenarios
    # (NB: some may be changed to 1 again after all -> fCNmax criterion)
    stage[is.na(BE)] <- 2
    
    # Define temporary variables for BE and stage for the second stage
    BE2 <- rep.int(NA, length(pes_tmp))
    s2 <- rep.int(2, length(pes_tmp))
    
    # Derive components for z-statistics
    Z11 <- qnorm(1 - p11) 
    Z12 <- qnorm(1 - p12)
    
    if (ssr.conditional == "no") {
      alpha_ssr <- cl$siglev[2]
      pwr_ssr <- targetpower
      lGMR_ssr <- if (usePE) pes_tmp else lGMR
    } else {
      # Derive conditional error rates
      alpha_ssr <- 1 - pnorm(pmin(
        (cl$cval[2] - sqrt(weight[1])*cbind(Z11, Z12)) / sqrt(1 - weight[1]),
        (cl$cval[2] - sqrt(weight[lw])*cbind(Z11, Z12)) / sqrt(1 - weight[lw])
      ))
      
      # Define target power for ssr
      pwr_ssr <- targetpower
      if ((ssr.conditional == "error_power") && (fCpower <= targetpower)) {
        # Use conditional estimated target power
        pwr_ssr <- 1 - (1 - targetpower) / (1 - pwr_s1)
      }
      
      if (usePE) {
        lGMR_ssr <- pes_tmp
      } else {
        # Set sign of lGMR to the sign of estimated point estimate
        # (Maurer et al call this 'adaptive planning step')
        sgn_pes_tmp <- ifelse(pes_tmp >= 0, 1, -1)
        lGMR_ssr <- abs(lGMR) * sgn_pes_tmp
      }
    }
    # Sample size for stage 2
    ptms <- proc.time()
    n2 <- .sampleN3(alpha = alpha_ssr, targetpower = pwr_ssr, 
                    ltheta0 = lGMR_ssr, mse = mses_tmp, 
                    ltheta1 = ltheta1, ltheta2 = ltheta2, method = pmethod)
    if (ssr.conditional == "no")
      n2 <- n2 - n1
    n2 <- pmax.int(pmin.int(n2, max.n - n1), min.n2)
    
    if (details) {
      cat("Time consumed (secs):\n")
      print(round((proc.time() - ptms), 1))
    }
    
    # Futility check regarding maximum overall sample size
    # - if n2 is infinite, then we consider this as fail too
    BE2[is.infinite(n2)] <- FALSE
    # - Otherwise check if n1+n2 is greater than fCNmax (if so, set to fail)
    if (fCrit != "no" && nms_match[3])
      BE2[n1 + n2 > fCNmax] <- FALSE
    s2[BE2 == FALSE] <- 1  # such a case is considered to be gone up to stage 1
    fut <- fut + if (all(is.na(BE2))) 0 else sum(BE2 == FALSE)
    
    # Carry over results from BE2 and s2 to BE, stage and ntot
    stage[is.na(BE)] <- s2
    BE[is.na(BE)] <- BE2
    n2 <- n2[is.na(BE2)]
    ntot[stage == 2] <- n1 + n2
    
    # Reduce relevant characteristics to the ones where stage 2 will be done
    p11 <- p11[is.na(BE2)]
    p12 <- p12[is.na(BE2)]
    Z11 <- Z11[is.na(BE2)]
    Z12 <- Z12[is.na(BE2)]
    
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
    
    rm(se.fac, df, n2, pes_s2, mses_s2, t1, t2)
   
    ## Combined evaluation via inverse normal approach
    # Derive further components for z-statistics 
    Z21 <- qnorm(1 - p21)
    Z22 <- qnorm(1 - p22)
    # For H01, test statistic at the end of second stage:
    Z01 <- pmax.int(
      sqrt(weight[1]) * Z11 + sqrt(1 - weight[1]) * Z21,
      sqrt(weight[lw]) * Z11 + sqrt(1 - weight[lw]) * Z21
    )
    # For H02, test statistic at the end of second stage:
    Z02 <- pmax.int(
      sqrt(weight[1]) * Z12 + sqrt(1 - weight[1]) * Z22,
      sqrt(weight[lw]) * Z12 + sqrt(1 - weight[lw]) * Z22
    )
    
    # Fill the remaining NAs with either TRUE or FALSE
    BE[is.na(BE)] <- (Z01 > cl$cval[2] & Z02 > cl$cval[2])
    # Could also do
    # BE[is.na(BE)] <- (p21 < alpha_ssr[, 1] & p22 < alpha_ssr[, 2])
    # i.e. test statistics are not really needed to be calculated.
    # Caution: p21 and p22 should however not be confused with the overall
    #          p-value resulting from this two-stage setting
    
    rm(Z11, Z21, Z12, Z22, Z01, Z02)
  }  # end if length(pes_tmp) > 0
  
  ### Define final output ------------------------------------------------------
  res <- list(
    # Information
    design = "2x2 crossover", method = "IN", alpha = cl$siglev, weight = weight,
    cval = cl$cval, max.comb.test = max.comb.test, n1 = n1, CV = CV, GMR = GMR, 
    usePE = usePE, targetpower = targetpower, fCpower = fCpower, 
    min.n2 = min.n2, max.n = max.n, ssr.conditional = ssr.conditional, 
    theta0 = theta0, theta1 = theta1, theta2 = theta2, 
    fCrit = fCrit, fCrange = c(fClower, fCupper), fCNmax = fCNmax, 
    pmethod = pmethod, nsims = nsims,
    # Results
    pBE = sum(BE) / nsims,
    pBE_s1 = sum(BE[ntot == n1]) / nsims,
    pct_stop_s1 = 100 * sum(ntot == n1) / nsims,
    pct_stop_fut = 100 * fut / nsims,
    pct_s2 = 100 * sum(ntot > n1) / nsims,
    nmean = mean(ntot),
    nrange = range(ntot),
    nperc = quantile(ntot, probs = npct)
  )
  
  if (details) {
    cat("Total time consumed (secs):\n")
    print(round((proc.time() - ptm), 1))
    cat("\n")
  }
  
  class(res) <- c("pwrtsd", "list")
  res
}  # end function power.2stage.in

# alias of the function
power.tsd.in <- power.2stage.in