### S3 print method for class evaltsd ------------------------------------------

print.evaltsd <- function(x, ...) {
  ### General information ------------------------------------------------------
  cat("TSD with 2x2 crossover\n")
  cat("Inverse Normal approach\n")
  if (x$max.comb.test) cat(" - maximum") else cat(" - standard")
  cat(" combination test\n")
  if (x$max.comb.test) {
    cat(" - weights (s1/s2) =", round(x$weight, 5), "\n")
  } else {
    cat(" - weight (s1) =", round(x$weight[[1]], 5), "\n")
  }
  cat(" - significance levels (s1/s2) =", round(x$alpha, 5), "\n")
  cat(" - critical values (s1/s2) =", round(x$cval, 5), "\n")
  cat(" - BE acceptance range = ", x$theta1," ... ", x$theta2, "\n", sep = "")

  ### Derived values -----------------------------------------------------------
  if (x$stage == 1) {
    not_pe <- if (x$usePE) "" else "not "
    cat(" - Observed point estimate from stage 1 is ", not_pe, "used for SSR\n",
        sep = "")
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
    cat("- Derived key statistics:\n")
    cat(sprintf("  z1 = %.5f, z2 = %.5f", x$z1, x$z2), ",\n", sep = "")
    cat("  Repeated CI = ",
        sprintf("(%.4f, %.4f)", x$RCI[[1]], x$RCI[[2]]), "\n", sep = "")

    if (x$stop_fut) {
      cat("- Futility criterion met:\n")
      if (x$futility[[1]] == 1) {
        cat("  * BE not declared and Power of first stage ",
            sprintf("(%.4f)", x$`Power Stage 1`), " > ",
            sprintf("%.4f", x$fCpower), "\n", sep = "")
      }
      if (x$futility[[2]] == 1) {
        fCrit <- x$fCrit
        if ("ci" %in% fCrit) {
          fCrit <- paste("90% CI", sprintf("(%.4f, %.4f)", x$CI90[[1]],
                                           x$CI90[[2]]))
        } else if ("pe" %in% fCrit) {
          fCrit <- "PE"
        }
        cat("  * ", fCrit," outside ", sprintf("%.4f", x$fCrange[1L]), " ... ", sprintf("%.4f", x$fCrange[2L]), "\n", sep = "")
      }
      if (x$futility[[3]] == 1) {
        cat("  * n2 = ", x$n2, " such that n1 + n2 > ", x$fCNmax, "\n", sep = "")
      }
      not_be <- if (x$stop_BE) "" else "not "
      cat("- Test for BE ", not_be,
          "positive (not considering any futility rule)\n", sep = "")
      cat("- Calculated n2 = ", x$n2, "\n", sep = "")
      cat("- Decision: Recommend to stop due to futility.\n")
    } else {
      cat("- No futility criterion met\n")
      not_be <- if (x$stop_BE) "" else "not "
      cat("- Test for BE ", not_be,
          "positive (not considering any futility rule)\n", sep = "")
      if (!x$stop_BE)
        cat("- Calculated n2 = ", x$n2, "\n", sep = "")
      cat("- Decision: ")
      if (x$stop_BE) {
        cat("Stop due to BE\n")
      } else {
        cat("Continue to stage 2 with ", x$n2, " subjects\n", sep = "")
      }
    }
  } else {
    cat("\nFinal analysis of second stage\n")
    cat("- Derived key statistics:\n")
    cat(sprintf("  z1 = %.5f, z2 = %.5f", x$z1, x$z2), ",\n", sep = "")
    cat("  Repeated CI = ",
        sprintf("(%.4f, %.4f)", x$RCI[[1]], x$RCI[[2]]), "\n", sep = "")
    cat("  Median unbiased estimate = ", sprintf("%.4f", x$MEUE), "\n",
        sep = "")
    cat("- Decision: ")
    if (x$stop_BE) {
      cat("BE achieved\n")
    } else {
      cat("BE not achieved\n")
    }
  }
}
