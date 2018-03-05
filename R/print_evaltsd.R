### S3 print method for class evaltsd ------------------------------------------

print.evaltsd <- function(x, ...) {
  ### General information ------------------------------------------------------
  cat("TSD with 2x2 crossover\n")
  cat("Inverse Normal approach\n")
  if (x$max.comb.test) cat(" - maximum") else cat(" - standard")
  cat(" combination test (")
  if (x$max.comb.test) cat("weights = ", x$weight[1], " ", x$weight[2], ")\n", 
                           sep = "") 
  else cat("weight = ", x$weight, ")\n", sep = "")
  cat(" - significance levels (s1/s2) =", round(x$cl$siglev, 5), "\n")
  cat(" - critical values (s1/s2) =", round(x$cl$cval, 5), "\n")
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
  ### Derived values -----------------------------------------------------------
  if (x$stage == 1) {
    cat("Interim analysis of first stage\n")
    ## TO DO
    ## Include t11, p11 etc.
    
    ## Futility regarding PE, CI, Nmax
    if ("no" %in% x$fCrit) {
      cat("No futility criterion regarding CI, PE or Nmax\n")
    } else {
      if (is.finite(x$fCNmax)) {
        cat("Futility criterion Nmax = ", x$fCNmax, ":", sep = "")
        if (x$futility[3] == 1) cat("met\n") else cat("not met\n")
      } else {
        cat("No futility criterion regarding Nmax\n")
      }
      if (x$fCrange[1L] > 0 && is.finite(x$fCrange[2L])) {
        fCrit <- x$fCrit
        fCrit <- if ("ci" %in% fCrit) "90% CI" else if ("pe" %in% fCrit) "PE"
        cat("Futility criterion ", fCrit," outside ", round(x$fCrange[1L], 4), 
            " ... ", round(x$fCrange[2L], 4), ":", sep = "")
        if (x$futility[2] == 1) {
          cat("met")
          if (!is.null(x$CI90$lower))
            cat("(90% CI = (", round(x$CI90$lower, 4), ", ", 
                round(x$CI90$upper, 4), "))", sep = "")
          cat("\n")
        } else { 
          cat("not met\n")
        }
      } else {
        cat("No futility criterion regarding PE or CI\n")
      }
    }
    ## Futility regarding Power of stage 1
    #cat("Power of stage 1 = ", round(x$`Power Stage 1`, 4), sep = "")
    #if (x$futility[1] == 1) cat(" > ", x$fCpower, "-> BE failed")
  } else {
    cat("Final analysis of second stage\n")
  }
}