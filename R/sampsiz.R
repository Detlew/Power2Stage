# -------------------------------------------------------------------------
# sample size function without all the overhead
# for 2x2 crossover or 2-group parallel
# power via nct approximation or exact
#
# author D. Labes
# -------------------------------------------------------------------------
# source("./R/sampsiz_n0.R")

.sampleN <- function(alpha=0.05, targetpower=0.8, ltheta0, ltheta1=log(0.8), 
                     ltheta2=log(1.25), mse, bk=2, method="nct")
{
  # return 'Inf' if ltheta0 not between or very near to ltheta1, ltheta2
  if ((ltheta0-ltheta1)<1.25e-5 | (ltheta2-ltheta0)<1.25e-5) {
    # debug
    # cat ("Inf returned for",exp(ltheta0),"\n")
    return(Inf)
  }

  # design characteristics for 2-group parallel and 2x2 crossover design
  steps <- 2     # stepsize for sample size search
  nmin  <- 4     # minimum n
  
  # log transformation assumed
  se    <- sqrt(mse)
  diffm <- ltheta0
  
  # start value from large sample approx. (hidden func.)
  # Jan 2015 changed to pure Zhang's formula
  # gives at least for 2x2 the best estimate (max diff to n: +-4)
  n  <- .sampleN0_3(alpha, targetpower, ltheta1, ltheta2, diffm, se, steps)
  if (n<nmin) n <- nmin
  if(method=="ls") return(n)
  df <- n-2 # both for 2x2 and parallel group
  pow  <- .calc.power(alpha, ltheta1, ltheta2, diffm, sem=se*sqrt(2/n), df, 
                      method)
  
  iter <- 0; imax <- 50
  # iter>50 is emergency brake
  # this is eventually not necessary, depends on quality of sampleN0
  # in experimentation I have seen max of 2-3 steps
  # reformulation with only one loop does not shorten the code considerable
  # --- loop until power <= target power, step-down
  down <- FALSE; up <- FALSE
  while (pow>targetpower) {
    if (n<=nmin) { 
      break
    }
    down <- TRUE
    n    <- n-steps     # step down if start power is to high
    iter <- iter+1
    df   <- n-2
    pow  <- .calc.power(alpha, ltheta1, ltheta2, diffm, sem=se*sqrt(2/n), df, 
                        method)
    
    if (iter>imax) break  
    # loop results in n with power too low
    # must step one up again. is done in the next loop
  }
  # --- loop until power >= target power
  while (pow<targetpower) {
    up   <- TRUE; down <- FALSE
    n    <- n+steps
    iter <- iter+1
    df   <- n-2
    pow  <- .calc.power(alpha, ltheta1, ltheta2, diffm, sem=se*sqrt(2/n), df, 
                        method)
    if (iter>imax) break 
  }
  nlast <- n
  if ((up & pow<targetpower) | (down & pow>targetpower) ) {
    n <- NA
  }
  
  #  return only n here
  return(n)
  
} # end of function
