# --------------------------------------------------------------------------
# power (or alpha) of 2-stage studies according to Potvin et al.
# methods "B" with 2 PK metrics
#
# Author D.L.
# --------------------------------------------------------------------------

power.tsd.2m <- function(alpha=c(0.0294,0.0294), n1, GMR, CV, targetpower=0.8,
                         pmethod=c("nct","exact", "shifted"), theta0, theta1, 
                         theta2, npct=c(0.05, 0.5, 0.95), nsims, setseed=TRUE, 
                         details=FALSE)
{
  if (missing(CV)) stop("CV's must be given.")
  if (any(CV<=0))  stop("CV's must be >0.")

  if (missing(n1)) stop("Number of subjects in stage 1 must be given.")
  if (n1<=0)       stop("Number of subjects in stage 1 must be >0.")

  if (length(alpha) != 2) stop("alpha must have two elements")

  if (missing(GMR)) GMR <- rep(0.95, 2)

  if (missing(theta1) & missing(theta2))  theta1 <- 0.8
  if (!missing(theta1) & missing(theta2)) theta2 <- 1/theta1
  if (missing(theta1) & !missing(theta2)) theta1 <- 1/theta2

  if (GMR[1]<=theta1 | GMR[1]>=theta2) stop("GMR[1] must be within acceptance range.")
  if (GMR[2]<=theta1 | GMR[2]>=theta2) stop("GMR[2] must be within acceptance range.")
  
  if (missing(theta0)) theta0 <- GMR
  else {
    if(length(theta0)!=2) stop("theta0 must have two elements")
  }

  if(missing(nsims)){
    nsims <- 1E5
    if(theta0[1] <= theta1 | theta0[1] >= theta2) nsims <- 1E6
    if(theta0[2] <= theta1 | theta0[2] >= theta2) nsims <- 1E6
  }

  # check if power calculation method is nct or exact
  pmethod <- match.arg(pmethod)

  if(details){
    cat(nsims,"sims. Stage 1")
  }
  # start timer
  ptm  <- proc.time()

  if (setseed) set.seed(1234567)

  ltheta1 <- log(theta1)
  ltheta2 <- log(theta2)
  lGMR    <- log(GMR)
  mlog    <- log(theta0)
  mse     <- CV2mse(CV)
  bk      <- 2   # 2x2x2 crossover design const
  # reserve memory
  BE      <- rep.int(NA, times=nsims)

# ----- stage 1 ----------------------------------------------------------
  Cfact <- bk/n1
  df    <- n1-2
  tval  <- qt(1-alpha[1], df)
  sdm   <- sqrt(mse*Cfact)
  # simulate point est. via normal distribution
  pes_m1   <- rnorm(n=nsims, mean=mlog[1], sd=sdm[1]) # metric 1, f.i. AUC
  pes_m2   <- rnorm(n=nsims, mean=mlog[2], sd=sdm[2]) # metric 2, f.i. Cmax
  # simulate mse via chi-squared distribution
  mses_m1  <- mse[1]*rchisq(n=nsims, df=df)/df
  mses_m2  <- mse[2]*rchisq(n=nsims, df=df)/df
  
  BE <- function(mses, pes)
  {
    hw <- tval*sqrt(Cfact*mses)
    loCL <- pes - hw
    upCL <- pes + hw
    BE <- loCL >= ltheta1 & upCL <= ltheta2
  }
  # make BE decision for stage 1
  BE_m1 <- BE(mses_m1, pes_m1)
  BE_m2 <- BE(mses_m2, pes_m2)
  
  # overall BE in stage 1
  BE <- BE_m1 & BE_m2
  
  # time for stage 1
  if(details){
    cat(" - Time consumed (secs):\n")
    print(round((proc.time()-ptm),1))
  }
  
  # BE == TRUE is decided yet
  # BE == FALSE not yet
  ntot <- rep(n1, nsims)
  if(sum(!BE)>0){
    ind <- !BE
    pes_m1  <- pes_m1[ind]
    mses_m1 <- mses_m1[ind]
    pes_m2  <- pes_m2[ind]
    mses_m2 <- mses_m2[ind]
    rm(ind)
    # 'abbreviated' method B, i.e. implicite power calculation via sample size 
    # estimation only
    # ------sample size for stage 2 -----------------------------------------
    # metric 1
    nt_m1 <- .sampleN2(alpha=alpha[2], targetpower=targetpower, ltheta0=lGMR[1],
                       mse=mses_m1, ltheta1=ltheta1, ltheta2=ltheta2,
                       method=pmethod)
    n2_m1 <- ifelse(nt_m1>n1, nt_m1 - n1, 0)
    # same for metric 2
    nt_m2 <- .sampleN2(alpha=alpha[2], targetpower=targetpower, ltheta0=lGMR[2],
                       mse=mses_m2, ltheta1=ltheta1, ltheta2=ltheta2,
                       method=pmethod)
    n2_m2 <- ifelse(nt_m2>n1, nt_m2 - n1, 0)
    
    n2 <- pmax(n2_m1, n2_m2 )
    
    # ------stage 2 evaluation ----------------------------------------------
    # simulate stage 2 and evalute the combined data from stage 1 + 2
    BE2 <- function(pes1, mses1, n2, nu)
    {
      # nu is number of metric
      #browser()
      m1    <- pes1
      SS1   <- (n1-2)*mses1
      nsim2 <- length(pes1)
      # to avoid warnings for n2=0 in rnorm() and rchisq()
      ow    <- options("warn")
      options(warn=-1)
      m2    <- ifelse(n2>0, rnorm(n=nsim2, mean=mlog[nu], sd=sqrt(mse[nu]*bk/n2)), 0)
      # ??? (n2-2) cancels out!
      SS2   <- ifelse(n2>2, (n2-2)*mse[nu]*rchisq(n=nsim2, df=n2-2)/(n2-2), 0)
      # reset options
      options(ow)
      SSmean <- ifelse(n2>0, (m1-m2)^2/(2/n1+2/n2), 0)
      nt     <- n1+n2
      df2    <- ifelse(n2>0, nt-3, n1-2)
      pe2    <- ifelse(n2>0, (n1*m1+n2*m2)/nt, pes1)
      mse2   <- ifelse(n2>0, (SS1+SSmean+SS2)/df2, mses1)
      # take care of memory
      rm(m1, m2, SS1, SS2, SSmean)
      # calculate CI for stage 2 with alpha[2]
      tval2 <- qt(1-alpha[2], df2)
      hw    <- tval2*sqrt(mse2*bk/nt)
      lower <- pe2 - hw
      upper <- pe2 + hw
      BE2   <- lower>=ltheta1 & upper<=ltheta2
      BE2    
    }
    # browser()
    BE2_m1 <- BE2(pes_m1, mses_m1, n2, nu=1)
    BE2_m2 <- BE2(pes_m2, mses_m2, n2, nu=2)
    # combine stage 1 & stage 2
    ntot[!BE] <- n2 + n1
    BE[!BE] <- BE2_m1 & BE2_m2
  }
    
  # the return list
  res <- list(design="2x2 crossover",
              method="B", alpha=alpha, VV=CV, n1=n1, GMR=GMR, 
              targetpower=targetpower, pmethod=pmethod,
              theta0=exp(mlog), theta1=theta1, theta2=theta2, nsims=nsims,
              # results
              pBE=sum(BE)/nsims,
              # Dec 2014 changed to
              pBE_s1=sum(BE[ntot==n1])/nsims,
              # was
              #pBE_s1=sum(BE[stage==2])/nsims,
              # Dec 2014 changed the meaning of pct_s2, was
              #pct_s2=100*length(BE[stage==2])/nsims,
              # whereby in case of unsymmetric alpha's stage 2 was also assigned
              # if n2=0 but not BE using alpha[1] to fascilate BE test with alpha[2]
              # now it is
              pct_s2=100*sum(ntot>n1)/nsims,
              # which simply means all those with n2>0
              nmean=mean(ntot), nrange=range(ntot), nperc=quantile(ntot, p=npct)
              #, ntot=ntot # experimental: return also all sample sizes
              )


  if (details){
    cat("Total time consumed (secs):\n")
    print(round((proc.time()-ptm),1))
    cat("\n")
  }

#  class(res) <- c("pwrtsd", "list")
  return(res)

} #end function
