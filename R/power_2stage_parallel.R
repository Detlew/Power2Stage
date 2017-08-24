# -----------------------------------------------------------------------------
# power (or alpha) of 2-stage studies with 2-group parallel design according to 
# Potvin et al. methods "B" and "C", modified to include a futility criterion Nmax,
# modified to use PE of stage 1 in sample size estimation
#
# author D.L.
# -----------------------------------------------------------------------------

power.2stage.p <- function(method=c("B","C"), alpha0=0.05, alpha=c(0.0294,0.0294), 
                           n1, GMR, CV, targetpower=0.8, 
                           pmethod=c("nct", "exact", "shifted"), 
                           usePE=FALSE, Nmax=Inf, test=c("welch", "t-test", "anova"),
                           theta0, theta1, theta2, npct=c(0.05, 0.5, 0.95),  
                           nsims, setseed=TRUE, details=FALSE)
{
  if (missing(CV)) stop("CV(s) must be given!")
  if (any(CV<=0))  stop("CV(s) must be >0!")
  # equal CV's
  if (length(CV)==1) {
    CVT <- CVR <-CV
  } else {
    CV <- CV[1:2] # truncate if longer then 2
    CVT <- CV[1]; CVR <- CV[2]
  }
  varT <- CV2mse(CVT)
  varR <- CV2mse(CVR)
  
  if (missing(n1))  stop("Number of subjects in stage 1 must be given!")
  if (length(n1)>2) stop("Length of n1 has to be 1 (scalar) or 2")
  if (any(n1<=0))   stop("Number of subjects in stage 1 must be >0!")
  if (length(n1)==1) {
    # total number given
    if (n1%%2!=0)  warning("Number of subjects in stage 1 should be even.\n",
                           "  Will be truncated to even.", immediate. = TRUE)
    n1 <- 2*trunc(n1/2)
    n1T <- n1R <- n1/2
  } else {
    n1T <- n1[1]
    n1R <- n1[2]
  }
  n1 <- n1T+n1R
  
  if (missing(GMR)) GMR <- 0.95
  
  if (missing(theta1) & missing(theta2))  theta1 <- 0.8
  if (!missing(theta1) & missing(theta2)) theta2 <- 1/theta1
  if (missing(theta1) & !missing(theta2)) theta1 <- 1/theta2
  
  if (GMR<=theta1 | GMR>=theta2) stop("GMR must be within acceptance range!")
  
  if (missing(theta0)) theta0 <- GMR
  
  if (n1>Nmax) stop("sum(n1)>Nmax doesn't make sense!")
  
  if(missing(nsims)){
    nsims <- 1E5
    if(theta0<=theta1 | theta0>=theta2) nsims <- 1E6
  }
  
  # check if Potvin B or C
  method  <- match.arg(method)
  # check test
  test  <- match.arg(test)
  # check if power calculation method is nct or exact
  pmethod <- match.arg(pmethod)
  
  if(details){
    cat(nsims,"sims. Stage 1")
  }
  # start timer
  ptm  <- proc.time()
  
  if (setseed) set.seed(1234567)

  # log transform
  ltheta1 <- log(theta1)
  ltheta2 <- log(theta2)
  lGMR    <- log(GMR)
  mlog    <- log(theta0)
  bk      <- 4   # 2-group parallel design constant
  bkni    <- 1   # design constant in terms of n(i)
  steps   <- 2   # for sample size search
  # reserve memory for the BE result
  BE      <- rep.int(NA, times=nsims)
  
# ----- stage 1 ----------------------------------------------------------
  
  nT <- n1T; nR <- n1R # short hand variables
  # do we need this?
  Cfact <- bk/n1
  df    <- n1T + n1R -2
  tval  <- qt(1-alpha[1], df) # ANOVA and t-test
  # simulate means (log domain) via normal distributions
  # astonishing enough the m1T are not the same, except the sign, 
  # if mlog=log(0.8) or mlog=log(1.25)
  # this is different from the situation of crossover where we simulate only
  # one mean, namely the difference
  m1T  <- rnorm(n=nsims, mean=mlog, sd=sqrt(varT/n1T))
  m1R  <- rnorm(n=nsims, mean=0, sd=sqrt(varR/n1R))
  # point est.
  pes   <- m1T-m1R
  # simulate variances via chi-squared distribution
  varsT  <- varT*rchisq(n=nsims, df=n1T-1)/(n1T-1)
  varsR  <- varR*rchisq(n=nsims, df=n1R-1)/(n1R-1)

  Vpooled <- ((n1T-1)*varsT + (n1R-1)*varsR)/(n1-2)
  se.fact <- sqrt(bkni*(1/n1T+1/n1R))

  if(method=="C"){
    # if method=C then calculate power for alpha0=0.05 and plan GMR
    # is this also the power of Welch t-test?
    # clear answer: no
    # see http://www.utdallas.edu/~ammann/stat6V99/node2.html
    # calculate CIs at alpha0 for all
    if (test=="t-test" || test=="anova"){
      pwr <- .calc.power(alpha=alpha0, ltheta1=ltheta1, ltheta2=ltheta2, 
                         diffm=lGMR, sem=sqrt(Vpooled)*se.fact, df=df,  
                         method=pmethod)
      hw    <- qt(1-alpha0, df)*sqrt(Vpooled)*se.fact
    } else {
      # Welch's t-test
      se  <- sqrt(varsT/nT + varsR/nR)
      dfs <- (varsT/nT + varsR/nR)^2/(varsT^2/nT^2/(nT-1)+varsR^2/nR^2/(nR-1))
      # here we use the Welch's df and se in power
      # necessary to truncate dfs? in case of pmethod=="exact"?
      # since OwensQOwen needs integer Df's
      if(pmethod=="exact") dfs <- trunc(dfs)
      pwr <- .calc.power(alpha=alpha0, ltheta1=ltheta1, ltheta2=ltheta2, 
                         diffm=lGMR, sem=se, df=dfs, method=pmethod)
      hw  <- qt(1-alpha0,dfs)*se
    }
    lower <- pes - hw
    upper <- pes + hw
    # fail or pass
    BE    <- lower>=ltheta1 & upper<=ltheta2
    # if power>0.8 then calculate CI for alpha=0.05
    # i.e. if power<0.8 then 
    BE[pwr<targetpower] <- NA # not yet decided
  }
  # method "B" or power<=0.8 in method "C"
  # then evaluate BE with alpha[1]
  Vpooled_tmp <- Vpooled[is.na(BE)]
  pes_tmp     <- pes[is.na(BE)]
  varsT_tmp   <- varsT[is.na(BE)] 
  varsR_tmp   <- varsR[is.na(BE)] 
  BE1 <- rep.int(NA, times=length(Vpooled_tmp))
  # calculate CI for alpha=alpha1
  if (test=="t-test" || test=="anova"){
    hw  <- tval*sqrt(Vpooled_tmp)*se.fact
  } else {
    # Welch's t-test
    se  <- sqrt(varsT_tmp/nT + varsR_tmp/nR)
    dfs <- (varsT_tmp/nT + varsR_tmp/nR)^2/(varsT_tmp^2/nT^2/(nT-1) + 
                                            varsR_tmp^2/nR^2/(nR-1))
    hw  <- qt(1-alpha[1], dfs)*se
    # we need dfs and se in following power monitoring step, so don't rm()!
  }
  rm(varsT_tmp, varsR_tmp)
  
  lower <- pes_tmp - hw
  upper <- pes_tmp + hw
  BE1   <- lower>=ltheta1 & upper<=ltheta2
  if (method=="C"){
    #if BE met -> PASS and stop
    #if not BE -> goto sample size estimation i.e flag BE1 as NA (not yet decided)
    BE1[!BE1] <- NA
  } else {
    # method B/E: evaluate power at alpha[2] and planGMR
    # should we have here also "B0"?
    if (test=="t-test" || test=="anova"){
      pwr <- .calc.power(alpha=alpha[2], ltheta1=ltheta1, ltheta2=ltheta2, 
                         diffm=lGMR, sem=sqrt(Vpooled_tmp)*se.fact, df=df, 
                         method=pmethod)
    } else {
      # here we use Welch's dfs and se from above
      # necessary to truncate dfs? in case of pmethod=="exact"?
      # since OwensQOwen needs integer Df's
      if(pmethod=="exact") dfs <- trunc(dfs)
      pwr <- .calc.power(alpha=alpha[2], ltheta1=ltheta1, ltheta2=ltheta2, 
                         diffm=lGMR, sem=se, df=dfs, method=pmethod)
    }  
    # this is decision scheme of MSDBE in case of alpha[1]!=alpha[2]
    # if BE met then decide BE regardless of power
    # if not BE and power<0.8 then goto stage 2
    #BE1[ !BE1 & pwr<targetpower] <- NA 
    
    # Xu et al., Potvin method E:
    # if not BE and if power >= 0.8 (targetpower) make a second BE evaluation 
    # with alpha[2]
    # only if alpha[1] != alpha[2] necessary, but works also without if(...)
    BE12  <- BE1 # reserve memory
    BE11  <- BE1
    # calculate CI for alpha=alpha2
    if (test=="t-test" || test=="anova"){
      hw <- qt(1-alpha[2], df)*sqrt(Vpooled_tmp)*se.fact
    } else {
      # Welch's t-test
      hw <- qt(1-alpha[2], dfs)*se
      # now we are done with them
      rm(dfs, se)
    }
    lower <- pes_tmp - hw
    upper <- pes_tmp + hw
    BE12 <- lower>=ltheta1 & upper<=ltheta2
    # if BE(a1) then BE1=TRUE, regardless of power
    BE1[BE11==TRUE] <- TRUE 
    # if not BE(a1) but power >= 0.8 then make BE decision at alpha2
    BE1[BE11==FALSE & pwr>=targetpower] <- BE12[BE11==FALSE & pwr>=targetpower]
    # if not BE(a1) and power<0.8 then not decided yet (marker NA)
    # will be further decided by futility criterion
    BE1[BE11==FALSE & pwr<targetpower]  <- NA 
    # keep care of memory
    rm(BE11, BE12)
  }
  # combine 'stage 0' from method C and stage 1
  BE[is.na(BE)] <- BE1
  # take care of memory
  # done with them
  rm(BE1, hw, lower, upper, pes_tmp, Vpooled_tmp)
  
  # time for stage 1
  if(details){
    cat(" - Time consumed (secs):\n")
    print(round((proc.time()-ptm),1))
  }

  # ------sample size for stage 2 -----------------------------------------
  ntot     <- rep(n1, times=nsims)
  stage    <- rep(1, times=nsims)
  # filter out those were stage 2 is necessary
  pes_tmp  <- pes[is.na(BE)]
  
  # Maybe we are already done with stage 1
  if(length(pes_tmp)>0){
    if(details){
      cat("Keep calm. Sample sizes for stage 2 (", length(pes_tmp),
          " studies)\n", sep="")
      cat("will be estimated. May need some time.\n")
    }
    # preliminary setting stage=2 for those not yet decided BE
    # may be altered for those with nts>Nmax or nts=Inf 
    # from sample size est. if pe outside acceptance range
    # see below
    stage[is.na(BE)] <- 2
    Vpooled_tmp <- Vpooled[is.na(BE)]
    m1T <- m1T[is.na(BE)]
    m1R <- m1R[is.na(BE)]
    varsT <- varsT[is.na(BE)]
    varsR <- varsR[is.na(BE)]
    BE2   <- rep.int(NA, times=length(Vpooled_tmp))
    s2    <- rep.int( 2, times=length(Vpooled_tmp))
    #------ sample size for stage 2 ---------------------------------------
    ptms <- proc.time()
    # sse always balanced
    # one correction step with Welch power in case of Welch's
    if (test=="anova") dfc <-"n-3" else dfc <- "n-2"
    if (usePE){
      # use mse1 & pe1 in sse like in the paper of Karalis/Macheras
      # sample size function returns Inf if pe1 is outside acceptance range
      # Aug. 2017: .sampleN2() uses now N-3 as df. was before N-2
      nts <- .sampleN2(alpha=alpha[2], targetpower=targetpower, ltheta0=pes_tmp,
                       mse=Vpooled_tmp, ltheta1=ltheta1, ltheta2=ltheta2, 
                       method=pmethod, bk=4, dfc=dfc)
      # in case of Welch's test the sample size may be too low
      # (se & sqrt(Vpooled)*se.fact are same if balanced, but df != dfs)
      # thus we raise them if necessary
      if (test=="welch"){
        # calculate power
        nT  <- nR <- nts/2
        se <- sqrt(varsT/nT + varsR/nR)
        dfs <- (varsT/nT + varsR/nR)^2/(varsT^2/nT^2/(nT-1) + varsR^2/nR^2/(nR-1))
        # in case of nts==Inf the dfs come out with NaN!
        # se is then zero and pwr comes out with NaN
        # we are setting the dfs to whatever not to thow errors
        dfs <- ifelse(!is.finite(nts), Inf, dfs)
        # necessary to truncate dfs? in case of pmethod=="exact"?
        # since OwensQOwen needs integer Df's
        if(pmethod=="exact") dfs <- trunc(dfs)
        pwr <- .calc.power(alpha=alpha[2], ltheta1=ltheta1, ltheta2=ltheta2, 
                           diffm=pes_tmp, sem=se, df=dfs, method=pmethod)
        nts <- ifelse(is.finite(nts) & pwr<targetpower, nts+steps, nts)
      }
    } else {
      # use mse1 & plan GMR to calculate sample size (original Potvin)
      nts <- .sampleN2(alpha=alpha[2], targetpower=targetpower, ltheta0=lGMR,
                       mse=Vpooled_tmp, ltheta1=ltheta1, ltheta2=ltheta2, 
                       method=pmethod, bk=4, dfc=dfc)
      # in case of Welch's test the sample size may be too low
      # thus we raise them if necessary
      if (test=="welch"){
        # calculate power for Welch's
        nT  <- nR <- nts/2
        se <- sqrt(varsT/nT + varsR/nR)
        dfs <- (varsT/nT + varsR/nR)^2/(varsT^2/nT^2/(nT-1) + varsR^2/nR^2/(nR-1))
        pwr <- .calc.power(alpha=alpha[2], ltheta1=ltheta1, ltheta2=ltheta2, 
                           diffm=lGMR, sem=se, df=dfs, method=pmethod)
        nts <- ifelse(is.finite(nts) & pwr<targetpower, nts+steps, nts)
      }
    }
    
    # The next is Jiri's strategy I think
    # nts are even, nts/2 for T and R (balanced at the end) 
    n2T <- ifelse(nts/2>n1T, nts/2 - n1T, 0)
    n2R <- ifelse(nts/2>n1R, nts/2 - n1R, 0)
    n2 <- n2T + n2R

    if(details){
      cat("Time consumed (secs):\n")
      print(round((proc.time()-ptms),1))
    }
    # futility rule: if nts > Nmax -> stay with stage 1 result not BE
    # ntotal = n1 reasonable?
    if (is.finite(Nmax) | any(!is.finite(nts))){
      # sample size may return Inf if PE is used in ss estimation
      # in that case we stay with stage 1
      BE2[!is.finite(n2) | (n1+n2)>Nmax] <- FALSE
      # and we are counting these for stage 1
      s2[BE2==FALSE]  <- 1
      # debug print
      # cat(sum(!BE2, na.rm=T)," cases with nts>Nmax or nts=Inf\n")
      # save 
      stage[is.na(BE)] <- s2
      # save the FALSE and NA in BE
      BE[is.na(BE)]    <- BE2
      # filter out those were BE was yet not decided
      pes_tmp  <- pes_tmp[is.na(BE2)]
      m1T <- m1T[is.na(BE2)]
      m1R <- m1R[is.na(BE2)]
      Vpooled_tmp <- Vpooled_tmp[is.na(BE2)]
      varsT <- varsT[is.na(BE2)]
      varsR <- varsR[is.na(BE2)]
      n2    <- n2[is.na(BE2)]
      n2T  <- n2T[is.na(BE2)]
      n2R  <- n2R[is.na(BE2)]
    } # end of futility Nmax

    # ----- simulate stage 2 data ------
    nsim2 <- length(pes_tmp)
    # to avoid warnings for ns2X==0 in rnorm() and ns2X<=1 in rchisq()
    ow    <- options("warn")
    options(warn=-1)
    m2T  <- ifelse(n2T>0, rnorm(n=nsim2, mean=mlog, sd=sqrt(varT/n2T)), 0)
    m2R  <- ifelse(n2R>0, rnorm(n=nsim2, mean=0, sd=sqrt(varR/n2R)), 0)
    # means T & R
    nT <- n1T+n2T
    nR <- n1R+n2R
    mT <- (n1T*m1T + n2T*m2T)/nT
    mR <- (n1R*m1R + n2R*m2R)/nR
    # point est.
    pes <- mT-mR
    # means for s1, s2 (stages)
    m1  <- (n1T*m1T + n1R*m1R)/(n1T+n1R)
    m2  <- ifelse((n2T+n2R)>0, (n2T*m2T + n2R*m2R)/(n2T+n2R),0)
    # total mean
    mt  <- (nT*mT + nR*mR)/(nT+nR)
    # simulate variances via chi-squared distribution
    # attention! in case of ns2X==1 rchisq gives NaN!
    # TODO: work out the correct way for ns2X==1
    # here we set it meanwhile to zero. result?
    vars2T  <- ifelse(n2T>1, varT*rchisq(n=nsim2, df=n2T-1)/(n2T-1), 0)
    vars2R  <- ifelse(n2R>1, varR*rchisq(n=nsim2, df=n2R-1)/(n2R-1), 0)
    # reset options
    options(ow) 
    
    # vars T/R over stage 1, stage 2
    # sy2 = sum of y squared
    # SQ = sum((y-mean)^2) = sum(y^2) - sum(y)^2/n = sum(y^2) - n*mean^2
    # var = SQ/(n-1)
    # bug in pre 0.4-5 (was mean^2/n instead of n*mean^2) which lead to
    # negativ residual variance and therefore BE = NA
    # example:
    # power.2stage.p(CV=0.4, n1=10, theta0=0.8, test="a")
    sy2T <- (n1T-1)*varsT + n1T*m1T^2
    sy2T <- ifelse(n2T>0, sy2T + (n2T-1)*vars2T + n2T*m2T^2, sy2T)
    sy2R <- (n1R-1)*varsR + n1R*m1R^2
    sy2R <- ifelse(n2R>0, sy2R + (n2R-1)*vars2R + n2R*m2R^2, sy2R)
    varsT <- (sy2T - nT*mT^2)/(nT-1)
    varsR <- (sy2R - nR*mR^2)/(nR-1)
    sy2   <- sy2T+sy2R
    sqtot <- (sy2-(nT+nR)*mt^2)
    rm(sy2T, sy2R, sy2)
    # calculate CI for stage 2 with alpha[2]
    if (test=="t-test" || test=="anova"){
      Vpooled <- ((nT-1)*varsT + (nR-1)*varsR)/(nT+nR-2)
      dfs <- (nT+nR-2)
      if (test=="anova"){
        # subtract stage effect
        # no subjects in stage 2 ((ns2R+ns2T)==0) may occure in case of 
        # high n1 and/or haybittle-peto alpha's
        # in the next expression Vpooled may come out as negative! why?
        # Example: CV=0.4, n1=10
        Vpooled <- ifelse(n2R+n2T>0,
                          (dfs*Vpooled - (n1*(m1-mt)^2 + n2*(m2-mt)^2))/(dfs-1),
                          Vpooled)
        # Vpooled according to potvin 2x2 crossover, Ben's Vermutung SSmean
        # but with (1/n1+1/n2) instead of (4/n1+4/n2)
        # Vp <- ifelse(n2R+n2T>0,
        #              (dfs*Vpooled - (m1-m2)^2/(1/n1+1/n2))/(dfs-1),
        #              Vpooled)
        # according to Anders paper
        # Vp <- ifelse(n2R+n2T>0,
        #              (sqtot - (n1*(m1-mt)^2 + n2*(m2-mt)^2)    # stage
        #               -(nT*(mT-mt)^2 + nR*(mR-mt)^2))/(dfs-1), # treatment
        #               Vpooled)
      }
      hw  <- qt(1-alpha[2],dfs)*sqrt(bkni*Vpooled*(1/nT+1/nR))
      rm(Vpooled, dfs)
    } else {
      # Welch's t-test
      se  <- sqrt(varsT/nT + varsR/nR)
      dfs <- (varsT/nT + varsR/nR)^2/(varsT^2/nT^2/(nT-1)+varsR^2/nR^2/(nR-1))
      hw  <- qt(1-alpha[2],dfs)*se
      rm(se, dfs)
    }
    # calculate CI and decide BE
    lower <- pes - hw
    upper <- pes + hw
    BE2   <- lower>=ltheta1 & upper<=ltheta2
    # combine stage 1 & stage 2
    ntot[is.na(BE)] <- n1+n2
    BE[is.na(BE)]   <- BE2
    # done with them
    rm(BE2, nts, lower, upper, hw, n2T, n2R, m2T, m2R)
  } # end stage 2 calculations
  # take care of memory
  rm(pes_tmp, pes)
  # the return list
  res <- list(design="2 parallel groups", method=method, 
              alpha0=ifelse(method=="C",alpha0,NA), alpha=alpha, CV=CV, 
              n1=c(n1T, n1R), GMR=GMR, test=test, targetpower=targetpower, 
              pmethod=pmethod, 
              theta0=exp(mlog), theta1=theta1, theta2=theta2, 
              usePE=usePE, Nmax=Nmax, nsims=nsims,
              # results 
              pBE=sum(BE)/nsims, pBE_s1=sum(BE[stage==1])/nsims,
              # Dec 2014: meaning of pct_s2 changed
              pct_s2=100*sum(ntot>n1)/nsims, 
              nmean=mean(ntot), nrange=range(ntot), 
              nperc=quantile(ntot, p=npct))
  
  # table object summarizing the discrete distri of ntot
  # only if usePE=FALSE or if usePE=TRUE then Nmax must be finite?
#  if (usePE==FALSE | (usePE==TRUE & is.finite(Nmax))){
    res$ntable <- table(ntot)
#  }
  
  if (details){
    cat("Total time consumed (secs):\n")
    print(round((proc.time()-ptm),1))
    cat("\n")
  }

  class(res) <- c("pwrtsd", "list")
  return(res)
  
} #end function
