up2even <- function(n) {
  if (n %% 2 != 0) return(2 * ceiling(n/2))
  n
}

# Check futility criterion for inverse normal methods --------------------------
futility_checks_in <- function(fCrit = "CI", fClower, fCupper, fCNmax, 
                               theta1, n1, min.n2 = 4) {
  stopifnot(is.character(fCrit))
  fCrit <- tolower(fCrit)
  fcrit_nms <- c("ci", "pe", "nmax", "no")  # correct possibilities
  nms_match <- fcrit_nms %in% fCrit  # check which fCrits are given
  if (sum(nms_match) == 0)
    stop("fCrit not correctly specified.")
  if (nms_match[4]) { # No futility criterion
    if (sum(nms_match[1:3]) > 0) {
      message("No futility will be applied.")
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
  list(fClower = fClower, fCupper = fCupper, fCNmax = fCNmax, 
       nms_match = nms_match)
}

# -----------------------------------------------------------------------------
# calculates the critical values for TSD evaluated via p-value combination
# Pocock method
critical.value.2stage <- function(alpha = 0.05, weight) {
  stopifnot(length(weight) <= 2, alpha > 0, alpha < 1)
  corr <- diag(1, nrow = 1 + length(weight), ncol = 1 + length(weight))
  if (length(weight) == 1) {
    varcov <- sqrt(weight)
  } else {
    r <- sqrt(prod(weight)) + sqrt(prod(1 - weight))
    varcov <- c(sqrt(weight), r)
  }
  corr[lower.tri(corr)] <- varcov
  corr[upper.tri(corr)] <- varcov
  
  f <- function(x) {
    # Could also use qmvnorm but is much slower
    pmvnorm(lower = rep.int(-Inf, nrow(corr)), upper = rep.int(x, nrow(corr)),
            mean = rep.int(0, nrow(corr)), corr = corr,
            algorithm = GenzBretz(maxpts = 100000, abseps = 1e-06)) - 
      (1 - alpha)
  }
  
  extint <- if (alpha > 0.5) "upX" else "no"
  z_crit <- uniroot(f, interval = c(0, qnorm(1 - alpha/(1 + length(weight)))), 
                    extendInt = extint)$root
  # critical value for stage 1 equal to stage 2 critical value
  z_crit <- c(z_crit, z_crit)
  list(cval = z_crit, siglev = 1 - pnorm(z_crit))
}

# -----------------------------------------------------------------------------
# utility function to create Wishart random draws of covariance matrices
# only dim(Sigma) == c(2,2) implemented
# attention! Sigma must be positiv definit!
# that prevents the use of rho=1 or rho=-1.
rWish2 <- function(n, df, Sigma)
{
  lendf <- length(df)
  # more checks to come
  stopifnot(dim(Sigma)==c(2,2))
  
  if (lendf==1) {
    ret <- rWishart(n, df, Sigma)
  } else {
    stopifnot(lendf==n)
    ret <- array(0, dim=c(2,2,n))
    for(i in seq_along(df)) {
      if(df[i] > 1) {
        ret[,,i] <- rWishart(1, df[i], Sigma)
      } else {
        # return null matrix
        # TODO: is this reasonable
        ret[,,i] <- matrix(0, nrow=2, ncol=2)
      }
    }
  }
  ret
}

# -----------------------------------------------------------------------------
# derive n, df and SEM (2x2x2 crossover setting)
get_n_df_sem <- function(n = NULL, df = NULL, mse, sem = NULL) {
  if (is.null(n)) {
    if (is.null(df)) {
      stop("Either n or df must be given.")
    } else {
      df <- df
      n <- df + 2 # only correct if no missing data and balanced
      sem <- if (is.null(sem)) sqrt(2 / n) * sqrt(mse) else sem
    }
  } else {
    if (length(n) == 1) {
      n <- if (is.finite(n)) PowerTOST:::nvec(n = n, grps = 2) else 
        rep(Inf, times = 2)
      if (n[1] != n[length(n)]) {
        message("Unbalanced design. n(i)=", paste(n, collapse="/"), " assumed.")
      } 
    } else {
      if (length(n) != 2) {
        stop("Length of n vector must be ", 2, "!")
      }
      if (any(n<1)) stop("All n(i) have to be >0.")
    }
    nc <- sum(1/n)
    n <- sum(n)
    se.fac <- sqrt(1/2 * nc)
    df <- if (is.null(df)) n - 2 else df
    sem <- if (is.null(sem)) se.fac * sqrt(mse) else sem
  }
  list(n = n, df = df, sem = sem)
}

#------------------------------------------------------------------------------
# function to check if function is called with .2stage. in its name
check2stage <- function(fname)
{
  if(grepl(".2stage.", fname)) 
    message("Function ", fname, " is deprecated.\n",
             "Use version with .tsd. in its name.")
}