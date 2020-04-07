# -------------------------------------------------------------------------
# sample size function without all the overhead
# for 2x2 crossover or 2-group parallel
# power via nct or shifted t approximation or exact
#
# Version with vectorize sample size w.r.t mse and/or ltheta0
#
# author D. Labes Jan 2015
# 
# Modified by BL, Jun 2017
# (vectorize wrt alpha)
# -------------------------------------------------------------------------
# is also used for parallel groups with bk=4

.sampleN3 <- function(alpha=0.05, targetpower=0.8, ltheta0, ltheta1=log(0.8), 
                      ltheta2=log(1.25), mse, method="nct", bk=2)
{

  # se and ltheta0/diffm must have the same length to vectorize propperly!
  if (length(mse)==1)     mse <- rep(mse, times=length(ltheta0))
  if (length(ltheta0)==1) ltheta0 <- rep(ltheta0, times=length(mse))
  if (length(targetpower)==1) targetpower <- rep(targetpower, times=length(mse))
  
  # Allow alpha to be scalar or a matrix.
  # If matrix, the interpretation is
  # - Require 2 columns
  # - Column 1 = alpha value for left hypothesis
  # - Column 2 = alpha value for right hypothesis
  # - Convention: If multiple rows, require same length as 
  #   max(length(diffm), length(sem)) and evaluate power element-wise 
  #   for each combination (alpha, diffm, sem)
  # If scalar, create a 1x2 matrix for consistency
  if (!is.matrix(alpha)) {
    if (length(alpha) == 1L) {
      alpha <- matrix(alpha, ncol = 2) # same alpha for both hypotheses
    } else {
      stop("alpha must be scalar (length 1) or matrix with 2 columns.")
    }
  }
  len <- max(length(ltheta0), length(mse))
  if (nrow(alpha) == 1) {
    # recycle
    alpha <- matrix(alpha, ncol = 2, nrow = len)
  } else {
    if (nrow(alpha) != len)
      stop("nrow(alpha) must be the same as length of delta1.")
  }
  dl <- ncol(alpha)
  if (dl > 2)
    stop("Number of columns of alpha should be 1 or 2.")
  
  # return 'Inf' if ltheta0 not between or very near to ltheta1, ltheta2
  #ns <- ifelse((ltheta0-ltheta1)<1.25e-5 | (ltheta2-ltheta0)<1.25e-5, Inf, 0)
  ns <- ifelse((ltheta0-ltheta1)<1e-04 | (ltheta2-ltheta0)<1e-04, Inf, 0)

  # design characteristics for 2-group parallel and 2x2 crossover design
  steps <- 2     # stepsize for sample size search
  nmin  <- 4     # minimum n
  
  se    <- sqrt(mse[is.finite(ns)])
  diffm <- ltheta0[is.finite(ns)]
  a <- alpha[is.finite(ns), , drop = FALSE]
  tp <- targetpower[is.finite(ns)]
  
  # start value from large sample approx. (hidden func.)
  # Jan 2015 changed to modified Zhang's formula
  # gives at least for 2x2 the best estimate (max diff to n: +-4)
  n <- .sampleN0_3(do.call(pmin, as.data.frame(a)), tp, ltheta1, ltheta2, 
                   diffm, se, steps, bk)
  n <- ifelse(n<nmin, nmin, n)
  
  # n may contain very large sample sizes which would result in very long
  # run time from loop below. Set all n's to Inf which are > 10^6
  ns[ns == 0] <- ifelse(n > 1e+06, Inf, 0)
  n <- n[n <= 1e+06]
  se <- sqrt(mse[is.finite(ns)])
  diffm <- ltheta0[is.finite(ns)]
  a <- alpha[is.finite(ns), , drop = FALSE]
  tp <- targetpower[is.finite(ns)]
  
  # method=="ls" is not used yet, in power.2stage.ssr() the 'original' ls approx
  # is used exactly as given in the paper of Golkowski et al.
  #if(method=="ls") return(n)
  
  # degrees of freedom n-2 for 2x2 crossover and 2-group parallel design 
  # since we have no stage term

  pow  <- .calc.power(a, ltheta1, ltheta2, diffm, sem=se*sqrt(bk/n), df=n-2, 
                      method)
  #iter <- rep(0, times=length(se)) 
  #imax <- 50
  index <- (pow > tp) & (n > nmin)
  while (any(index)) {
    n[index] <- n[index] - steps
    pow_tmp <- .calc.power(a[index, , drop=FALSE], ltheta1, ltheta2, 
                           diffm[index], sem=se[index]*sqrt(bk/n[index]), 
                           df=n[index]-2, method)
    pow[index] <- pow_tmp
    index <- (pow > tp) & (n > nmin)
  }
  index <- (pow < tp)
  while (any(index)) {
    n[index] <- n[index] + steps
    pow_tmp <- .calc.power(a[index, , drop=FALSE], ltheta1, ltheta2, 
                           diffm[index], sem=se[index]*sqrt(bk/n[index]), 
                           df=n[index]-2, method)
    pow[index] <- pow_tmp
    index <- (pow < tp)
  }
  
  # combine the Inf and n
  ns[ns==0] <- n
  n <- ns
  #  return only n here
  return(n)
  
} # end of function