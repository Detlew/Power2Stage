final.2stage.in <- function(alpha, weight, max.comb.test = TRUE, GMR1, CV1, n1,
                            df1 = NULL, SEM1 = NULL, GMR2, CV2, n2, df2 = NULL,
                            SEM2 = NULL, theta1, theta2) {

  # Check if called with .2stage. version
  check2stage(fname=as.character(sys.call())[1])

  ### Error handling and default value set-up ----------------------------------
  if (missing(GMR1) || missing(GMR2))
    stop("GMR1 and GMR2 must be given.")
  if (missing(CV1) || missing(CV2))
    stop("CV1 and CV2 must be given.")
  if (CV1 <= 0 || CV2 <= 0)
    stop("CV1 and CV2 must be > 0.")
  if (missing(n1) || missing(n2))
    stop("n1 and n2 must be given.")
  if (n1 < 3 || n2 < 3)
    stop("n1 and n2 must be at least 3.")
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

  if (missing(theta1) && missing(theta2))  theta1 <- 0.8
  if (!missing(theta1) && missing(theta2)) theta2 <- 1/theta1
  if (missing(theta1) && !missing(theta2)) theta1 <- 1/theta2
  lGMR1   <- log(GMR1)
  lGMR2   <- log(GMR2)
  ltheta1 <- log(theta1)
  ltheta2 <- log(theta2)
  mse1    <- CV2mse(CV1)
  mse2    <- CV2mse(CV2)
  des <- get_n_df_sem(n = n1, df = df1, mse = mse1, sem = SEM1)
  n1 <- des$n
  df1 <- des$df
  sem1 <- des$sem
  des <- get_n_df_sem(n = n2, df = df2, mse = mse2, sem = SEM2)
  n2 <- des$n
  df2 <- des$df
  sem2 <- des$sem

  ### Calculate adjusted critical levels ---------------------------------------
  cl <- if (length(alpha) == 1) critical.value.2stage(alpha, weight) else
    list(cval = qnorm(1 - alpha), siglev = alpha)

  ### Evaluation of Stage 2 ----------------------------------------------------
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
  ## Bioequivalence after stage 2?
  BE <- (Z01 > cl$cval[2] && Z02 > cl$cval[2])

  ## Calculate corresponding exact repeated CI
  rci <- repeated_ci(diff1 = lGMR1, diff2 = lGMR2, sem1 = sem1, sem2 = sem2,
                     df1 = df1, df2 = df2, a1 = cl$siglev[1], a2 = cl$siglev[2],
                     weight = weight, stage = 2)

  ## Calculate median unbiased estimate (Section 8.3.3 in Wassmer & Brannath)
  # Calculate conservative estimate: For an equivalence setting this means
  # to calculate lower and upper bound, and then take the worst case
  meue <- vector("numeric", 2)
  meue[[1]] <- median_unbiased_pe(diff1 = lGMR1, diff2 = lGMR2, sem1 = sem1,
                                  sem2 = sem2, df1 = df1, df2 = df2,
                                  a1 = cl$siglev[1], a0 = 1,
                                  weight = weight, lower_bnd = TRUE)

  meue[[2]] <- median_unbiased_pe(diff1 = lGMR1, diff2 = lGMR2, sem1 = sem1,
                                  sem2 = sem2, df1 = df1, df2 = df2,
                                  a1 = cl$siglev[1], a0 = 1,
                                  weight = weight, lower_bnd = FALSE)

  idx_max <- which.max(abs(meue))
  exp_meue <- exp(meue[[idx_max]])

  ### Define final output ------------------------------------------------------
  res <- list(
    stage = 2L, alpha = cl$siglev, cval = cl$cval, weight = weight,
    max.comb.test = max.comb.test, GMR1 = GMR1, CV1 = CV1, n1 = as.integer(n1),
    df1 = df1, SEM1 = sem1, GMR2 = GMR2, CV2 = CV2, n2 = as.integer(n2),
    df2 = df2, SEM2 = sem2, theta1 = theta1, theta2 = theta2,
    z1 = Z01, z2 = Z02, RCI = exp(rci), MEUE = exp_meue, stop_BE = BE
  )
  class(res) <- c("evaltsd", "list")
  res
}  # end of function final.2stage.in

# alias of function
final.tsd.in <- final.2stage.in
