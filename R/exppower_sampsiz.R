#------------------------------------------------------------------------------
# Author: dlabes
#------------------------------------------------------------------------------
# Approximate "expected" power according to Julious book
# taking into account the uncertainty of an estimated se with 
# dfse degrees of freedom
# Only for log-transformed data.
# Raw function: see the call in exppower.TOST()
.exppower <- function(alpha=0.05, ltheta1, ltheta2, diffm, sem, dfse, df)
{
  tval <- qt(1 - alpha, df, lower.tail = TRUE)
  d1   <- sqrt((diffm-ltheta1)^2/sem^2)
  d2   <- sqrt((diffm-ltheta2)^2/sem^2)
  # in case of diffm=ltheta1 or =ltheta2 and se=0
  # d1 or d2 have then value NaN (0/0)
  d1[is.nan(d1)] <- 0
  d2[is.nan(d2)] <- 0
  
  pow  <- pt(d1,dfse,tval) + pt(d2,dfse,tval) - 1
  
  return(ifelse(pow>=0,pow,0))
  
}  


# -------------------------------------------------------------------------
# sample size function without all the overhead
# for 2x2 crossover or 2-group parallel
# power via Julious expected power
#
# Version with vectorized sample size w.r.t mse and/or ltheta0
#
# author D. Labes Feb 2015
# -------------------------------------------------------------------------
# library(PowerTOST)
# source("./R/sampsiz_n0.R")
# is also used for parallel groups with bk=4

.expp.sampleN <- function(alpha=0.05, targetpower=0.8, ltheta0, ltheta1=log(0.8), 
                          ltheta2=log(1.25), mse, dfmse, bk=2)
{
  
  # se and ltheta0/diffm must have the same length to vectorize propperly!
  if (length(mse)==1)     mse <- rep(mse, times=length(ltheta0))
  if (length(ltheta0)==1) ltheta0 <- rep(ltheta0, times=length(mse))
  if (length(dfmse)==1)   dfmse <- rep(dfmse, times=length(mse))
  
  # return 'Inf' if ltheta0 not between or very near to ltheta1, ltheta2
  ns <- ifelse((ltheta0-ltheta1)<1.25e-5 | (ltheta2-ltheta0)<1.25e-5, Inf, 0)
  
  # design characteristics for 2-group parallel and 2x2 crossover design
  # df for the design as an unevaluated expression
  #dfe   <- parse(text="n-2", srcfile=NULL)
  steps <- 2     # stepsize for sample size search
  nmin  <- 4     # minimum n
  
  se     <- sqrt(mse[is.finite(ns)])
  diffm  <- ltheta0[is.finite(ns)]
  
  # start value from large sample approx. (hidden func.)
  # Jan 2015 changed to pure Zhang's formula
  # gives at least for 2x2 the best estimate (max diff to n: +-4)
  n <- .sampleN0_3(alpha, targetpower, ltheta1, ltheta2, diffm, se, steps, bk)
  n <- ifelse(n<nmin, nmin, n)
  
  df <- n-2
  pow <- .exppower(alpha=alpha, ltheta1=ltheta1, ltheta2=ltheta2, 
                   diffm=diffm, sem=se*sqrt(bk/n), dfse=dfmse, df=df)
  iter <- rep(0, times=length(se)) 
  imax <- 50
  # iter>50 is emergency brake
  # this is eventually not necessary, depends on quality of sampleN0
  # in experimentation I have seen max of 2-3 steps
  # reformulation with only one loop does not shorten the code considerable
  # --- loop until power <= target power, step-down
  #  down <- FALSE; up <- FALSE
  index <- pow>targetpower & n>nmin
  #browser()
  while (any(index)) {
    #    down <- TRUE
    n[index]    <- n[index]-steps     # step down if start power is to high
    iter[index] <- iter[index]+1
    df   <- n-2
    pow[index] <- .exppower(alpha=alpha, ltheta1=ltheta1, ltheta2=ltheta2, 
                            diffm=diffm[index], sem=se[index]*sqrt(bk/n[index]), 
                            dfse=dfmse[index], df=df[index])
    
    index <- pow>targetpower & n>nmin & iter<=imax
    # loop results in n with power too low
    # must step one up again. is done in the next loop
  }
  
  # --- loop until power >= target power
  index <- pow<targetpower
  while (any(index)) {
    #    up   <- TRUE; down <- FALSE
    n[index] <- n[index]+steps
    iter[index] <- iter[index]+1
    df   <- n-2
    pow[index] <- .exppower(alpha=alpha, ltheta1=ltheta1, ltheta2=ltheta2, 
                            diffm=diffm[index], sem=se[index]*sqrt(bk/n[index]), 
                            dfse=dfmse[index], df=df[index])
    index <- pow<targetpower & iter<=imax
  }
  #   nlast <- n
  #   if ((up & pow<targetpower) | (down & pow>targetpower) ) {
  #     n <- NA
  #   }
  # combine the Inf and n
  ns[ns==0] <- n
  n <- ns
  #  return only n here
  return(n)
  
} # end of function
