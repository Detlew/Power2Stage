### S3 print method for class evaltsd ------------------------------------------

print.evaltsd <- function(x, ...) {
  ### General information ------------------------------------------------------
  cat("TSD with 2x2 crossover\n")
  cat("Inverse Normal approach\n")
  if (x$max.comb.test) cat(" - maximum") else cat(" - standard")
  cat(" combination test\n")
  cat(" - significance levels (s1/s2) =", round(x$alpha, 5), "\n")
  cat(" - critical values (s1/s2) =", round(x$cval, 5), "\n")
  ### Derived values -----------------------------------------------------------
  if (x$stage == 1) {
    if (x$ssr.conditional == "no") {
      cat(" - without conditional error rates and conditional (estimated target) power\n")
    } else {
      cat(" - with ")
      if (x$ssr.conditional == "error") {
        cat("conditional error rates\n")
      } else {
        if (x$fCpower > x$targetpower)
          cat("conditional error rates\n")
        else
          cat("conditional error rates and conditional (estimated target) power\n")
      }
    }
    cat("\nInterim analysis of first stage\n")
    if (x$stop_s1) {
      if (x$stop_fut) {
        cat("- Stop due to futility:\n")
        if (x$futility[[1]] == 1) {
          cat("  BE not declared and Power of first stage (",
              round(x$`Power Stage 1`, 4), ") > ", round(x$fCpower, 4), "\n",
              sep = "")
        } else if (x$futility[[2]] == 1) {
          fCrit <- x$fCrit
          if ("ci" %in% fCrit) {
            fCrit <- paste("90% CI", sprintf("(%.4f, %.4f)", x$CI90$lower, x$CI90$upper))
          } else if ("pe" %in% fCrit) {
            fCrit <- "PE"
          }
          cat("  ", fCrit," outside ", round(x$fCrange[1L], 4),  " ... ",
              round(x$fCrange[2L], 4), "\n", sep = "")
        } else if (x$futility[[3]] == 1) {
          cat("  n2 such that n1 + n2 > ", x$fCNmax, "\n", sep = "")
        }
      } else {
        cat("- Stop because BE can be concluded\n")
        cat("- Derived key statistics:\n")
        cat(sprintf("  p11 = %.5f, p12 = %.5f", x$p11, x$p12), "\n")
        cat("  Repeated CI = ",
            sprintf("(%.4f, %.4f)", x$RCI$lower, x$RCI$upper), "\n", sep = "")
      }
    } else {
      cat("- BE not achieved\n")
      cat("- No futility criterion met\n")
      cat("- Continue to stage 2 with n2 = ", x$n2, " subjects\n", sep = "")
    }
  } else {
    cat("\nFinal analysis of second stage\n")
    cat("- BE acceptance range = ", x$theta1," ... ", x$theta2,"\n", sep="")
    if (x$BE) {
      cat("- BE achieved\n")
    } else {
      cat("- BE not achieved\n")
    }
    cat("- Derived key statistics:\n")
    cat(sprintf("  z01 = %.5f, z02 = %.5f", x$z01, x$z02), "\n")
    cat("  Repeated CI = ",
        sprintf("(%.4f, %.4f)", x$RCI$lower, x$RCI$upper), "\n", sep = "")
    cat("  Median unbiased estimate = ", sprintf("%.4f", x$MEUE), "\n", sep = "")
  }
}
