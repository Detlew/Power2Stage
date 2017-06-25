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
# library(PowerTOST)
# source("./R/sampsiz_n0.R")
# is also used for parallel groups with bk=4

.sampleN3 <- function(alpha=0.05, targetpower=0.8, ltheta0, ltheta1=log(0.8), 
                      ltheta2=log(1.25), mse, method="nct", bk=2)
{

  # se and ltheta0/diffm must have the same length to vectorize propperly!
  if (length(mse)==1)     mse <- rep(mse, times=length(ltheta0))
  if (length(ltheta0)==1) ltheta0 <- rep(ltheta0, times=length(mse))
  
  # Allow alpha to be scalar or a matrix. If matrix, the interpretation is
  # - Require 2 columns
  # - Column 1 = alpha value for left hypothesis
  # - Column 2 = alpha value for right hypothesis
  # - Convention: If multiple rows require diffm & sem to be of the same length 
  #   and evaluate power element-wise for each combination (alpha, diffm, sem)
  if (is.atomic(alpha) && !is.matrix(alpha)) {
    # We enforce matrix structure -> recycling will not work, so do it manually
    if (length(alpha) != max(length(ltheta0), length(mse))) {
      alpha <- rep.int(alpha, max(length(ltheta0), length(mse)))
    }
    alpha <- matrix(alpha, ncol = 1)
  } else { 
    if (nrow(alpha) != length(ltheta0) || length(ltheta0) != length(mse))
      stop("number of rows of alpha must match length of diffm and sem.")
  }
  dl <- ncol(alpha)
  if (dl > 2)
    stop("Number of columns of alpha should be 1 or 2.")
  
  # return 'Inf' if ltheta0 not between or very near to ltheta1, ltheta2
  ns <- ifelse((ltheta0-ltheta1)<1.25e-5 | (ltheta2-ltheta0)<1.25e-5, Inf, 0)

  # design characteristics for 2-group parallel and 2x2 crossover design
  steps <- 2     # stepsize for sample size search
  nmin  <- 4     # minimum n
  
  se    <- sqrt(mse[is.finite(ns)])
  diffm <- ltheta0[is.finite(ns)]
  
  # start value from large sample approx. (hidden func.)
  # Jan 2015 changed to modified Zhang's formula
  # gives at least for 2x2 the best estimate (max diff to n: +-4)
  n <- .sampleN0_3(do.call(pmin, as.data.frame(alpha)),
                   targetpower, ltheta1, ltheta2,  diffm, se, steps, bk)
  n <- ifelse(n<nmin, nmin, n)
  
  # method=="ls" is not used yet, in power.2stage.ssr() the 'original' ls approx
  # is used exactly as given in the paper of Golkowski et al.
  if(method=="ls") return(n)
  
  # degrees of freedom as expression
  # n-2 for 2x2 crossover and 2-group parallel design
  ##dfe <- parse(text="n-2", srcfile=NULL)  # is that needed?
  # or should that read n-3? see Kieser/Rauch
  #dfe <- parse(text="n-3", srcfile=NULL)
  
  #df   <- eval(dfe)
  pow  <- .calc.power(alpha, ltheta1, ltheta2, diffm, sem=se*sqrt(bk/n), df=n-2, 
                      method)
  #iter <- rep(0, times=length(se)) 
  #imax <- 50
  index <- (pow > targetpower) & (n > nmin)
  while (any(index)) {
    n[index] <- n[index] - steps
    pow_tmp <- .calc.power(alpha[index, , drop=FALSE], ltheta1, ltheta2, 
                           diffm[index], sem=se[index]*sqrt(bk/n[index]), 
                           df=n[index]-2, method)
    pow[index] <- pow_tmp
    index <- (pow > targetpower) & (n > nmin)
  }
  index <- (pow < targetpower)
  while (any(index)) {
    n[index] <- n[index] + steps
    pow_tmp <- .calc.power(alpha[index, , drop=FALSE], ltheta1, ltheta2, 
                           diffm[index], sem=se[index]*sqrt(bk/n[index]), 
                           df=n[index]-2, method)
    pow[index] <- pow_tmp
    index <- (pow < targetpower)
  }
  
  # combine the Inf and n
  ns[ns==0] <- n
  n <- ns
  #  return only n here
  return(n)
  
} # end of function