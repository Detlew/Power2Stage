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