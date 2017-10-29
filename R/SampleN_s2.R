# ------------------------------------------------------------------------------
# sample size for stage 2 for the two-stage design scheme
# based on the inverse normal method (Maurer et al.)
#
# author D. Labes based on Ben's power.2stage.in()
# ------------------------------------------------------------------------------

SampleNs2.in <- function(alpha, weight, max.comb.test = TRUE, n1, CV1, GMR1,
                         targetpower = 0.8, theta0, theta1, theta2, 
                         usePE = FALSE, min.n2 = 4, 
                         ssr.conditional = c("error_power", "error", "no"),
                         pmethod = c("nct", "exact", "shifted") )
{
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
  if (missing(CV1))  stop("CV1 must be given!")
  if (CV1 <= 0)      stop("CV1 must be >0!")
  if (targetpower <= 0 || targetpower >= 1) 
    stop("targetpower must be within (0, 1)")
  if (min.n2 < 4) 
    stop("min.n2 has to be at least 4.")
  if (min.n2 %% 2 != 0) {  # make it even
    min.n2 <- 2 * floor(min.n2 / 2) + 2
    message("min.n2 rounded up to next even integer", min.n2)
  }
  if (missing(GMR1)) stop("GMR1 must be given.")
  if (missing(theta0)) theta0 <- 0.95
  if (missing(theta1) && missing(theta2))  theta1 <- 0.8
  if (!missing(theta1) && missing(theta2)) theta2 <- 1/theta1
  if (missing(theta1) && !missing(theta2)) theta1 <- 1/theta2
  if (GMR1 <= theta1 || GMR1 >= theta2) 
    stop("GMR1 must be within acceptance range!")
  if (theta0 <= theta1 || theta0 >= theta2) 
    stop("theta0 must be within acceptance range!")
  
  # Check usage of conditional error rates / power
  ssr.conditional <- match.arg(ssr.conditional)
  
  # Check if power calculation method is nct, exact or shifted
  pmethod <- match.arg(pmethod)

  ltheta0 <- log(theta0)
  ltheta1 <- log(theta1)
  ltheta2 <- log(theta2)
  lGMR    <- log(GMR1)
  mse     <- CV2mse(CV1)
  bk      <- 2  # 2x2x2 crossover design const in terms of N (total)
  
  ### Calculate adjusted critical levels ---------------------------------------
  # length(alpha)=1 => if length(weight)=2, max.comb.test=TRUE is assumed here
  cl <- if (length(alpha) == 1) critical.value.2stage(alpha, weight) else
    # if two alphas, ignore possibly given weights and only use those alphas
    list(cval = qnorm(1 - alpha), siglev = alpha)

  se.fac <- sqrt(bk / n1)
  df <- n1 - 2

  # The 2 one-sided test statistics
  t1 <- (lGMR - ltheta1) / sqrt(mse) / se.fac
  t2 <- (lGMR - ltheta2) / sqrt(mse) / se.fac
  
  # p-values of first stage
  p11 <- pt(t1, df = df, lower.tail = FALSE)
  p12 <- pt(t2, df = df, lower.tail = TRUE)
  
  # TODO: return n2=0 if p11 & p12 significant
  
  # power of stage 1
  pwr_s1 <- .calc.power(alpha = cl$siglev[1], ltheta1 = ltheta1, 
                        ltheta2 = ltheta2, diffm = lGMR, 
                        sem = se.fac * sqrt(mse), df = df, method = pmethod)
  # TODO: return n2=0 pwr_s1 > targetpower. 
  # or better with fCpower argument which is not there at present?

  # Derive components for z-statistics
  Z11 <- qnorm(1 - p11) 
  Z12 <- qnorm(1 - p12)
  
  if (ssr.conditional == "no") {
    alpha_ssr <- cl$siglev[2]
    pwr_ssr <- targetpower
    lGMR_ssr <- if (usePE) lGMR else ltheta0
  } else {
    # Derive conditional error rates
    alpha_ssr <- 1 - pnorm(pmin(
      (cl$cval[2] - sqrt(weight[1])*cbind(Z11, Z12)) / sqrt(1 - weight[1]),
      (cl$cval[2] - sqrt(weight[lw])*cbind(Z11, Z12)) / sqrt(1 - weight[lw])
    ))
    
    # Define target power for ssr
    pwr_ssr <- targetpower
    if (ssr.conditional == "error_power") {
      # Use conditional power
      # this formula may give negative values if pwr_s1 >= targetpower. 
      # TODO: what then?
      pwr_ssr <- 1 - (1 - targetpower) / (1 - pwr_s1)
      if (pwr_ssr<0) pwr_ssr <- 0 # ??? with this setting n2=Inf returns
    }
    if (usePE) {
      lGMR_ssr <- lGMR
    } else {
      # Set sign of lGMR to the sign of estimated point estimate
      # (Maurer et al call this 'adaptive planning step')
      sgn_pe <- ifelse(lGMR >= 0, 1, -1)
      lGMR_ssr <- abs(ltheta0) * sgn_pe
    }
  }
  # Sample size for stage 2
  n2 <- .sampleN3(alpha = alpha_ssr, targetpower = pwr_ssr, 
                  ltheta0 = lGMR_ssr, mse = mse, ltheta1 = ltheta1, 
                  ltheta2 = ltheta2, method = pmethod)
  if (ssr.conditional == "no") n2 <- n2 - n1
  
  n2 <- pmax.int(n2, min.n2)
  
  browser()
  
  # what should we return? only n2?
  # or something (data.frame) like sampleN.TOST() or sampleN2.TOST() returns?
  n2
  
}
  