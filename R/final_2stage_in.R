final.2stage.in <- function(GMR1, GMR2, CV1, CV2, n1, n2, df1 = NULL, df2 = NULL,
                            SEM1 = NULL, SEM2 = NULL, alpha, weight,
                            max.comb.test = TRUE, theta1, theta2) {
  ### TO DO: Error handling
  if (missing(alpha))
    alpha <- 0.05
  if (length(alpha) > 2) 
    stop("Length of alpha must be <= 2.")
  if (!missing(weight) && length(weight) > 2)
    stop("Length of weight must be <= 2.")
  if (length(alpha) == 1 && missing(weight))
    weight <- if (max.comb.test) c(0.5, 0.25) else 0.5
  if (length(alpha) == 2) {
    if (missing(weight)) {
      stop("weight must be specified.")
    } else {
      message(paste0("Note: Adjusted alphas are specified, it is assumed that",
                     " the specified weight(s) are in line with the alpha",
                     " values."))
    }
  }
  lw <- length(weight)
  if (max.comb.test) {
    if (lw != 2)
      stop("Two weights, w and w*, are required for maximum combination test.")
  } else {
    if (lw != 1)
      stop("One weight, w, is required for standard combination test.")
  }
  if (any(weight <= 0) || any(weight >= 1))
    stop("Weight(s) not properly specified, must be > 0 and < 1.")
  if (any(alpha <= 0) || any(alpha >= 1))
    stop("Alpha(s) not properly specified, must be > 0 and < 1.")
  
  lGMR1   <- log(GMR1)
  lGMR2   <- log(GMR2)
  ltheta1 <- log(theta1)
  ltheta2 <- log(theta2)
  mse1    <- CV2mse(CV1)
  mse2    <- CV2mse(CV2)
  cl <- if (length(alpha) == 1) critical.value.2stage(alpha, weight) else
    list(cval = qnorm(1 - alpha), siglev = alpha)
  
  if (missing(n1)) {
    if (is.null(df1)) {
      stop("Either n1 or df1 must be given.")
    } else {
      df1 <- df1
      n1 <- df1 + 2
      sem1 <- if (is.null(SEM1)) sqrt(2 / n1) * sqrt(mse1) else SEM1
    }
  } else {
    df1 <- if (is.null(df1)) n1 - 2 else df1
    sem1 <- if (is.null(SEM1)) sqrt(2 / n1) * sqrt(mse1) else SEM1
  }
  if (missing(n2)) {
    if (is.null(df2)) {
      stop("Either n2 or df2 must be given.")
    } else {
      df2 <- df2
      n2 <- df2 + 2
      sem2 <- if (is.null(SEM2)) sqrt(2 / n2) * sqrt(mse2) else SEM2
    }
  } else {
    df2 <- if (is.null(df2)) n2 - 2 else df2
    sem2 <- if (is.null(SEM2)) sqrt(2 / n2) * sqrt(mse2) else SEM2
  }
  
  t11 <- (lGMR1 - ltheta1) / sem1
  t12 <- (lGMR1 - ltheta2) / sem1
  p11 <- pt(t11, df = df1, lower.tail = FALSE)
  p12 <- pt(t12, df = df1, lower.tail = TRUE)
  
  t21 <- (lGMR2 - ltheta1) / sem2
  t22 <- (lGMR2 - ltheta2) / sem2
  p21 <- pt(t21, df = df2, lower.tail = FALSE)
  p22 <- pt(t22, df = df2, lower.tail = TRUE)
  
  Z11 <- qnorm(1 - p11) 
  Z12 <- qnorm(1 - p12)
  Z21 <- qnorm(1 - p21)
  Z22 <- qnorm(1 - p22)
  
  # For H01, test statistic at the end of second stage:
  Z01 <- pmax.int(
    sqrt(weight[1]) * Z11 + sqrt(1 - weight[1]) * Z21,
    sqrt(weight[lw]) * Z11 + sqrt(1 - weight[lw]) * Z21
  )
  # For H02, test statistic at the end of second stage:
  Z02 <- pmax.int(
    sqrt(weight[1]) * Z12 + sqrt(1 - weight[1]) * Z22,
    sqrt(weight[lw]) * Z12 + sqrt(1 - weight[lw]) * Z22
  )
  
  BE <- (Z01 > cl$cval[2] & Z02 > cl$cval[2])
  
  ## Calculate corresponding exact repeated CI
  f <- function(t) {
    Z11 <- qnorm(1 - pt((lGMR1 - t) / sem1, df = df1, lower.tail = FALSE))
    Z21 <- qnorm(1 - pt((lGMR2 - t) / sem2, df = df2, lower.tail = FALSE))
    pmax.int(
      sqrt(weight[1]) * Z11 + sqrt(1 - weight[1]) * Z21,
      sqrt(weight[lw]) * Z11 + sqrt(1 - weight[lw]) * Z21
    ) - cl$cval[2]
  }
  g <- function(t) { 
    Z12 <- qnorm(1 - pt((lGMR1 - t) / sem1, df = df1, lower.tail = TRUE))
    Z22 <- qnorm(1 - pt((lGMR2 - t) / sem2, df = df2, lower.tail = TRUE))
    pmax.int(
      sqrt(weight[1]) * Z12 + sqrt(1 - weight[1]) * Z22,
      sqrt(weight[lw]) * Z12 + sqrt(1 - weight[lw]) * Z22
    ) - cl$cval[2]
  }
  search_int <- lGMR2 + c(-5, 5) * sem2
  ll <- uniroot(f, interval = search_int)$root
  lu <- ll + 2 * (lGMR2 - ll)  # corresponding upper bound
  ru <- uniroot(g, interval = search_int)$root
  rl <- ru - 2 * (ru - lGMR1)  # corresponding lower bound
  # CI for equivalence problem via intersection
  lower <- max(ll, rl)
  upper <- min(lu, ru)
  ci <- if (upper < lower) NA else list(lower = exp(lower), upper = exp(upper))
  
  # TO DO: expand
  res <- list(
    eRCI = ci, BE_s2 = BE
  )
  res
}  # end of function final.2stage.in

# alias of function
final.tsd.in <- final.2stage.in