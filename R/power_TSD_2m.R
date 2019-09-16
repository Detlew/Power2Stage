# --------------------------------------------------------------------------
# power (or alpha) of 2-stage studies according to Potvin et al. methods "B"
# with 2 PK metrics
#
# Author D.L.
# --------------------------------------------------------------------------

power.tsd.2m <- function(alpha=c(0.0294,0.0294), CV, n1, rho=0, GMR,
                         targetpower=0.8, pmethod=c("nct","exact", "shifted"),
                         powerstep=FALSE, theta0, theta1, theta2,
                         npct=c(0.05, 0.5, 0.95), nsims, setseed=TRUE, details=FALSE)
{
  if(missing(CV))  stop("CVs must be given.")
    else {
      if(length(CV)==1) CV <- rep(CV,2)
      if(length(CV)!=2) stop("GMR must have 2 elements.")
      if(any(CV<=0))   stop("CVs must be >0.")
    }
  if(missing(n1)) stop("Number of subjects in stage 1 must be given.")
  if(length(n1)!=1) {
    n1 <- n1[1]
    message("n1 must be scalar. First element of n1 used.")
    if(n1<=0) stop("Number of subjects in stage 1 must be >0.")
  }
  if(length(alpha)==1) {
    alpha <- rep(alpha, 2)
    message("Scalar alpha used for both stages.")
  }
  if(length(alpha) != 2) stop("alpha must have two elements.")

  if(missing(GMR)) GMR <- rep(0.95, 2)
    else {
      if(length(GMR)==1) GMR <- rep(GMR, 2)
      if(length(GMR)!=2) stop("GMR must have 2 elements.")
    }

  if (missing(theta1) & missing(theta2))  theta1 <- 0.8
  if (!missing(theta1) & missing(theta2)) theta2 <- 1/theta1
  if (missing(theta1) & !missing(theta2)) theta1 <- 1/theta2
  stopifnot(length(theta1)==1, length(theta2)==1)

  if(any(GMR<=theta1) | any(GMR>=theta2)) stop("GMRs must be within acceptance range.")

  if (missing(theta0)) theta0 <- GMR
  else {
    if(length(theta0)==1) theta0 <- rep(theta0, 2)
    if(length(theta0)!=2) stop("theta0 must have two elements.")
  }

  #if(rho!=0) warning("rho != 0 is only experimental.", call. = FALSE)
  stopifnot(length(rho)==1, rho >= -1, rho <= 1)

  if (rho < -1 || rho > 1)
    stop("Correlation must be within {-1, +1}.")

  # prevent error in rWishart(n, df, Sigma)
  if (rho == -1) rho <- -1 + .Machine $double.eps
  if (rho == +1) rho <- +1 - .Machine $double.eps

  if(missing(nsims)){
    nsims <- 1E5
    if(any(theta0<=theta1) | any(theta0>=theta2)) nsims <- 1E6
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
  pes  <- matrix(0, nrow=nsims, ncol=2)
  mses <- matrix(0, nrow=nsims, ncol=2)
  if (rho==0){
    pes[ , 1] <- rnorm(n=nsims, mean=mlog[1], sd=sdm[1]) # metric 1, f.i. AUC
    pes[ , 2] <- rnorm(n=nsims, mean=mlog[2], sd=sdm[2]) # metric 2, f.i. Cmax
    # simulate mse via chi-squared distribution
    mses[ , 1] <- mse[1]*rchisq(n=nsims, df=df)/df
    mses[ , 2] <- mse[2]*rchisq(n=nsims, df=df)/df
  } else {
    # multivariate normal with rho
    # TODO: check this
    #browser()
    # (population) variance-covariance matrix of means
    s_m <- diag(sdm^2)
    s_m[1,2] <- s_m[2,1] <- rho*s_m[1,1]*s_m[2,2]
    pes <- rmvnorm(nsims, mean=mlog, sigma=s_m)
    # simulate covariance matrices via Wishart distribution
    s_mse  <- diag(mse)
    s_mse[1,2] <- s_mse[2,1] <- rho*sqrt(s_mse[1,1]*s_mse[2,2])
    #or do we have here df=n-1?
    covm      <- rWish2(n=nsims, df=df, Sigma=s_mse)/df
    mses[, 1] <- covm[1, 1, ]
    mses[, 2] <- covm[2, 2, ]
    # take care of memory
    rm(covm)
  }

  BE <- function(nu)
  {
    hw <- tval*sqrt(Cfact*mses[, nu])
    loCL <- pes[, nu] - hw
    upCL <- pes[, nu] + hw
    BE <- loCL >= ltheta1 & upCL <= ltheta2
    BE
  }
  # make BE decision for stage 1
  BE_m1 <- BE(nu=1)
  BE_m2 <- BE(nu=2)

  # overall BE in stage 1
  BE <- BE_m1 & BE_m2
  #browser()
  # set BE to NA if BE not decided yet
  BE[!BE] <- NA

  if(powerstep) {
    # metric 1 power
    pwr_m1 <- .calc.power(alpha=alpha[2], ltheta1=ltheta1, ltheta2=ltheta2,
                          diffm=lGMR[1], sem=sqrt(bk*mses[, 1]/n1), df=df,
                          method=pmethod)
    BE[BE_m2 & !BE_m1 & pwr_m1>=targetpower] <- FALSE
    # metric 2 power
    pwr_m2 <- .calc.power(alpha=alpha[2], ltheta1=ltheta1, ltheta2=ltheta2,
                          diffm=lGMR[2], sem=sqrt(bk*mses[, 2]/n1), df=df,
                          method=pmethod)
    BE[BE_m1 & !BE_m2 & pwr_m2>=targetpower] <- FALSE
    # take care of memory
    rm(pwr_m1, pwr_m2)
  }
  # take care of memory
  rm(BE_m1, BE_m2)
  # time for stage 1
  if(details){
    cat(" - Time consumed (secs):\n")
    print(round((proc.time()-ptm),1))
  }

  # BE == TRUE/FALSE is decided yet
  # BE == NA not yet decided
  ntot <- rep(n1, nsims)
  if(sum(is.na(BE)) > 0) {
    ptms <- proc.time()
    ind  <- is.na(BE)
    pes  <- pes[ind, ]
    mses <- mses[ind, ]
    rm(ind)

    # ------sample size for stage 2 -----------------------------------------
    if(GMR[1]==GMR[2]){
      if(details){
        cat("Keep calm. Sample sizes for stage 2 (", sum(is.na(BE)),
            " studies)\n", sep="")
        cat("will be estimated. May need some time.\n")
      }
      # use the max. mse for sample size
      mses_max <- pmax(mses[, 1], mses[, 2])
      nt_max <- .sampleN2(alpha=alpha[2], targetpower=targetpower, ltheta0=lGMR[1],
                          mse=mses_max, ltheta1=ltheta1, ltheta2=ltheta2,
                          method=pmethod)
      n2 <- ifelse(nt_max>n1, nt_max - n1, 0)
      # keep care of memory
      rm(mses_max, nt_max)
    } else {
      if(details){
        cat("Keep calm. Sample sizes for stage 2 (2x ", sum(is.na(BE)),
            " studies)\n", sep="")
        cat("will be estimated. May need some time.\n")
      }

      # metric 1
      nt_m1 <- .sampleN2(alpha=alpha[2], targetpower=targetpower, ltheta0=lGMR[1],
                         mse=mses[, 1], ltheta1=ltheta1, ltheta2=ltheta2,
                         method=pmethod)
      n2_m1 <- ifelse(nt_m1>n1, nt_m1 - n1, 0)
      # same for metric 2
      nt_m2 <- .sampleN2(alpha=alpha[2], targetpower=targetpower, ltheta0=lGMR[2],
                         mse=mses[, 2], ltheta1=ltheta1, ltheta2=ltheta2,
                         method=pmethod)
      n2_m2 <- ifelse(nt_m2>n1, nt_m2 - n1, 0)
      n2 <- pmax(n2_m1, n2_m2 )
      rm(nt_m1, n2_m1, nt_m2, n2_m2)
    }
    if(details){
      cat("Time consumed (secs):\n")
      print(round((proc.time()-ptms),1))
    }
    #browser()
    # ------stage 2 evaluation ----------------------------------------------
    # simulate stage 2 and evalute the combined data from stage 1 + 2
    nsim2 <- nrow(pes)
    pes2  <- matrix(0, ncol=2, nrow=nsim2 )
    SS2   <- matrix(0, ncol=2, nrow=nsim2 )
    ow    <- options("warn")
    options(warn=-1)
    #browser()
    if (rho!=0){
      # we are simulating for mean=0, sigma = matrix 1 rho, rho 1
      sigma2 <- diag(1,2,2)
      sigma2[1,2] <- sigma2[2,1] <- rho
      pes2 <- rmvnorm(nsim2, mean=rep(0,2), sigma=sigma2)
      # now we transform to the actual variance and mean
      # TODO: is this correct? Or have we to include rho in some way
      # for the transformation?
      pes2[, 1] <- ifelse(n2>0, pes2[ ,1]*sqrt(mse[1]*bk/n2) + mlog[1], 0)
      pes2[, 2] <- ifelse(n2>0, pes2[ ,2]*sqrt(mse[2]*bk/n2) + mlog[2], 0)
      # SS2 via Wishart distribution
      #browser()
      covm     <- rWish2(n=nsim2, df=n2-2, Sigma=s_mse)
      SS2[, 1] <- covm[1, 1, ]
      SS2[, 2] <- covm[2, 2, ]
    } else {
      # pe's of stage 2 data via independent normal distri's
      pes2[,1] <- ifelse(n2>0, rnorm(n=nsim2, mean=mlog[1], sd=sqrt(mse[1]*bk/n2)), 0)
      pes2[,2] <- ifelse(n2>0, rnorm(n=nsim2, mean=mlog[2], sd=sqrt(mse[2]*bk/n2)), 0)
      # SS2 via independent chi-squared distribution
      SS2[, 1]   <- ifelse(n2>2, (n2-2)*mse[1]*rchisq(n=nsim2, df=n2-2)/(n2-2), 0)
      SS2[, 2]   <- ifelse(n2>2, (n2-2)*mse[2]*rchisq(n=nsim2, df=n2-2)/(n2-2), 0)
    }
    options(ow)

    BE2 <- function(nu)
    {
      # this function sees pes[], mses[], pes2[], SS2[] and n2[]
      # nu is number of metric
      #browser()
      m1     <- pes[, nu]
      SS1    <- (n1-2)*mses[, nu]
      nsim2  <- length(m1)
      m2     <- pes2[, nu]
      SSmean <- ifelse(n2>0, (m1-m2)^2/(2/n1+2/n2), 0)
      nt     <- n1+n2
      df2    <- ifelse(n2>0, nt-3, n1-2)
      pe2    <- ifelse(n2>0, (n1*m1+n2*m2)/nt, m1)
      mse2   <- ifelse(n2>0, (SS1+SSmean+SS2[, nu])/df2, mses[, nu])
      # take care of memory
      rm(m1, m2, SS1, SSmean)
      # calculate CI for stage 2 with alpha[2]
      tval2  <- qt(1-alpha[2], df2)
      hw     <- tval2*sqrt(mse2*bk/nt)
      lower  <- pe2 - hw
      upper  <- pe2 + hw
      BE2    <- lower>=ltheta1 & upper<=ltheta2
      BE2
    }
    # browser()
    BE2_m1 <- BE2(nu=1)
    BE2_m2 <- BE2(nu=2)
    # combine stage 1 & stage 2
    ntot[is.na(BE)] <- n2 + n1
    BE[is.na(BE)]   <- BE2_m1 & BE2_m2
  }

  # the return list
  res <- list(design="2x2 crossover",
              method="B2m", alpha=alpha, CV=CV, n1=n1, GMR=GMR, rho=rho,
              targetpower=targetpower, pmethod=pmethod,
              theta0=exp(mlog), theta1=theta1, theta2=theta2, usePE=FALSE,
              nsims=nsims,
              # results
              pBE=sum(BE)/nsims,
              pBE_s1=sum(BE[ntot==n1])/nsims,
              pct_s2=100*sum(ntot>n1)/nsims,
              # which simply means all those with n2>0
              nmean=mean(ntot), nrange=range(ntot), nperc=quantile(ntot, p=npct)
              #, ntot=ntot # experimental: return also all sample sizes
              )
  res$ntable <- table(ntot)

  if (details){
    cat("Total time consumed (secs):\n")
    print(round((proc.time()-ptm),1))
    cat("\n")
  }

  class(res) <- c("pwrtsd", "list")
  return(res)

} #end function
