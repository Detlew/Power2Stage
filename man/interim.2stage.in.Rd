\encoding{utf-8}
\name{interim.2stage.in}
\alias{interim.2stage.in}
\title{
Interim analysis of first stage data of 2-stage 2x2 crossover designs based on
the Inverse Normal method
}
\description{
  Following the design scheme according to \code{power.2stage.in} the function
  performs the interim analysis of the first stage data.
}
\usage{
interim.2stage.in(alpha, weight, max.comb.test = TRUE, targetpower = 0.8,
                  GMR1, n1, CV1, df1 = NULL, SEM1 = NULL, theta1, theta2, 
                  GMR, usePE = FALSE, min.n2 = 4, max.n = Inf, 
                  fCpower = targetpower, fCrit = "CI", fClower, fCupper, fCNmax, 
                  ssr.conditional = c("error_power", "error", "no"), 
                  pmethod = c("nct", "exact", "shifted"))
}
\arguments{
  \item{alpha}{
    If one element is given, the overall one-sided significance level (not
    the adjusted level for \ifelse{html}{\out{stage&nbsp;1}}{stage 1}). If two elements are given, the adjusted 
    one-sided alpha levels for \ifelse{html}{\out{stage&nbsp;1}}{stage 1} and
    \ifelse{html}{\out{stage&nbsp;2}}{stage 2}, respectively.\cr
    If missing, defaults to \verb{0.05}.
  }
  \item{weight}{
    Weight(s) of \ifelse{html}{\out{stage&nbsp;1}}{stage 1}. Note that using
    the notation from Maurer et al weight corresponds to information fraction,
    other literature may refer to sqrt(weight) as being the weight.
    \verb{weight} must either contain one element (in case of 
    \verb{max.comb.test = FALSE}) or two elements (in case of
    \verb{max.comb.test = TRUE}).\cr
    If missing, defaults to \verb{0.5} for \verb{max.comb.test = FALSE} and to
    \verb{c(0.5, 0.25)} for \verb{max.comb.test = TRUE}.
  }
  \item{max.comb.test}{
    Logical; if \verb{TRUE} (default) the maximum combination test will be
    used, otherwise the standard combination test.
  }
  \item{targetpower}{
    Desired (overall) target power to declare \acronym{BE} at the end of the trial.
  }
  \item{GMR1}{
    Observed ratio of geometric means (T/R) of \ifelse{html}{\out{stage&nbsp;1}}{stage 1} data.
  }
  \item{n1}{
    Sample size of \ifelse{html}{\out{stage&nbsp;1}}{stage 1}.
  }
  \item{CV1}{
    Observed coefficient of variation of the intra-subject variability of
    \ifelse{html}{\out{stage&nbsp;1}}{stage 1} (use \emph{e.g.}, 0.3 for 30\%).
  }
  \item{df1}{
    Optional; Error degrees of freedom of 
    \ifelse{html}{\out{stage&nbsp;1}}{stage 1} that can be specified instead
    of n1.
  }
  \item{SEM1}{
    Optional; Standard error of the difference of means of 
    \ifelse{html}{\out{stage&nbsp;1}}{stage 1} that can be specified instead
    of CV1. Must be on additive scale (i.e. usually log-scale).
  }
  \item{theta1}{
    Lower bioequivalence limit. Defaults to 0.8.
  }
  \item{theta2}{
    Upper bioequivalence limit. Defaults to 1.25.
  }
  \item{GMR}{
    Assumed ratio of geometric means (T/R) to be used in power calculation
    for \ifelse{html}{\out{stage&nbsp;1}}{stage 1} and sample size re-estimation
    for \ifelse{html}{\out{stage&nbsp;2}}{stage 2}.
  }
  \item{usePE}{
    If \verb{TRUE} the sample size re-estimation and the power of \ifelse{html}{\out{stage&nbsp;1}}{stage 1} 
    futility inspection is done with the observed point estimate 
    (\acronym{PE}) of the treatment difference in 
    \ifelse{html}{\out{stage&nbsp;1}}{stage 1}.\cr
    Defaults to \verb{FALSE}.
  }
  \item{min.n2}{
    Minimum sample size of \ifelse{html}{\out{stage&nbsp;2}}{stage 2}. Defaults to 4.\cr
    If the sample size re-estimation step gives a sample size for
    \ifelse{html}{\out{stage&nbsp;2}}{stage 2} less than \verb{min.n2}, then
    \verb{min.n2} will be used for \ifelse{html}{\out{stage&nbsp;2}}{stage 2}.
  }
  \item{max.n}{
    Maximum overall sample size \ifelse{html}{\out{stage&nbsp;1}}{stage 1} +
    \ifelse{html}{\out{stage&nbsp;2}}{stage 2}.\cr
    This is \emph{not} a futility criterion regarding the maximum sample size! If
    \verb{max.n} is set to a finite value and the sample size re-estimation gives a
    sample size for \ifelse{html}{\out{stage&nbsp;2}}{stage 2} (\verb{n2}) such
    that \verb{n1 + n2 > max.n}, then the sample size for
    \ifelse{html}{\out{stage&nbsp;2}}{stage 2} will be set to \verb{n2 = max.n - n1}.\cr
    Defaults to \verb{Inf}, \emph{i.e.}, no constraint on the re-estimated sample size.
  }
  \item{fCpower}{
    Threshold for power monitoring step to decide on futility for cases where
    \acronym{BE} has not been achieved after
    \ifelse{html}{\out{stage&nbsp;1}}{stage 1}: If \acronym{BE} has not been
    achieved after \ifelse{html}{\out{stage&nbsp;1}}{stage 1} and the power for
    \ifelse{html}{\out{stage&nbsp;1}}{stage~1} is greater than or equal to
    \verb{fCpower}, then the study will be considered a failure.\cr\cr
    See \sQuote{Details} for more information on the choice of
    \verb{fCpower}.
  }
  \item{fCrit}{
    Futility criterion to use: \verb{"No"} (no futility criterion regarding
    observed point estimate, confidence interval and maximum sample size), 
    \verb{"PE"} (observed point estimate of the geometric mean ratio from 
    \ifelse{html}{\out{stage&nbsp;1}}{stage~1}), 
    \verb{"CI"} (90\% confidence interval of the geometric mean ratio from
    \ifelse{html}{\out{stage&nbsp;1}}{stage~1}), \code{"Nmax"} (overall maximum
    sample size); or a combination thereof (concatenate abbreviations; see
    \sQuote{Examples}). Defaults to "CI".
  }
  \item{fClower}{
    Lower futility limit for the \acronym{PE} or \acronym{CI} of
    \ifelse{html}{\out{stage&nbsp;1}}{stage~1}.\cr
    If the \acronym{PE} or \acronym{CI} is completely outside of \verb{fClower}
    \ldots \verb{fCupper} the study is stopped due to futility (not \acronym{BE}).\cr
    May be missing. If \verb{"PE"} or \verb{"CI"} is specified within \verb{fCrit},
    the default will be set to 0.8 for \verb{fCrit = "PE"} or 0.95 for \verb{fCrit = "CI"}.
    If neither \verb{"PE"} nor \verb{"CI"} is specified within \verb{fCrit}, there
    will be no futility constraint regarding point estimate or confidence interval
    from \ifelse{html}{\out{stage&nbsp;1}}{stage~1} (regardless of any 
    specification of \verb{fClower} and/or \verb{fCupper}).
  }
  \item{fCupper}{
    Upper futility limit for the \acronym{PE} or \acronym{CI} of
    \ifelse{html}{\out{stage&nbsp;1}}{stage 1}.\cr
    Analogous to \verb{fClower}: Will be set to \verb{1/fClower} if missing.
  }
  \item{fCNmax}{
    Futility criterion regarding maximum sample size. If the determined sample size
    for \ifelse{html}{\out{stage&nbsp;2}}{stage 2} (\verb{n2}) is such that
    \verb{n1 + n2 > fCNmax}, the study will not continue to
    \ifelse{html}{\out{stage&nbsp;2}}{stage 2} and stopped due to futility (not
    \acronym{BE}).\cr
    If \verb{"Nmax"} is specified within \verb{fCrit} and argument \verb{fCNmax}
    is missing, the value will be set to \verb{fCNmax = 4*n1}. If \verb{"Nmax"} is
    not specified within \verb{fCrit}, then there will be no futility constraint
    regarding maximum sample size (regardless of any specification of \verb{fCNmax}).
  }
  \item{ssr.conditional}{
    Method for sample size re-estimation step: \verb{"no"} does not use
    conditional error rates nor the conditional estimated target power for 
    the second stage, \verb{"error"} uses conditional error rates for the 
    second stage, and \verb{"error_power"} uses both conditional error rates 
    and the conditonal estimated target power for the second stage.\cr
    Defaults to \verb{"error_power"}.\cr\cr
    See also \sQuote{Details}.
  }
  \item{pmethod}{
    Power calculation method, also to be used in the sample size estimation for
    \ifelse{html}{\out{stage&nbsp;2}}{stage 2}.\cr
    Implemented are \verb{"nct"} (approximate calculations via non-central
    \emph{t}-distribution, \verb{"exact"} (exact calculations via Owen\enc{’}{'}s Q),
    and \verb{"shifted"} (approximate calculation via shifted central \emph{t}-distribution
    like in the paper of Potvin \emph{et al.})\cr
    In contast to \code{power.2stage.in} the default value here is \verb{"exact"}. 
  }
}
\details{
  If \verb{ssr.conditional = "error_power"}, the design scheme generally 
  calculates the conditional estimated target power of the second stage and 
  uses this value as desired target power in the sample size re-estimation 
  process. If \code{fCpower} > \code{targetpower}, then the conditional power
  may actually be negative. This does not seem sensible. Therefore, for such
  cases the desired target power for the sample size re-calculation will be set 
  to \code{targetpower}.
}
\value{
  Returns an object of class \verb{"evaltsd"} with all the input arguments and results
  as components. As part of the input arguments a component \code{cval} is also
  presented, containing the critical values for \ifelse{html}{\out{stage&nbsp;1}}{stage 1} and 2 according to the
  input based on \code{alpha}, \code{weight} and \code{max.comb.test}.\cr
  The class \verb{"evaltsd"} has an S3 print method.\cr\cr
  The results are in the components:
  \item{futility}{Three dimensional vector with either 0 or 1. The first
  component represents futility due to Power of first stage > \code{fCpower},
  the second futility due to \verb{CI} (or \verb{PE}) outside of 
  \code{fClower ... fCupper}, the third futility due to \code{n1 + n2 > fCNmax}.}
  \item{t11}{Observed value of test statistic for first hypothesis.}
  \item{t12}{Observed value of test statistic for second hypothesis.}
  \item{p11}{Observed p-value for first hypothesis.}
  \item{p12}{Observed p-value for second hypothesis.}
  \item{CI90}{List with lower and upper 90\% CI limits for observed ratio
  of geometric means from \ifelse{html}{\out{stage&nbsp;1}}{stage 1}. If \verb{fCrit != "CI"} results will be \code{NULL}.}
  \item{Power Stage 1}{Calculated power of \ifelse{html}{\out{stage&nbsp;1}}{stage 1}.}
  \item{n2}{Required (total) sample size for \ifelse{html}{\out{stage&nbsp;2}}{stage 2} (zero if study does not
  continue to \ifelse{html}{\out{stage&nbsp;2}}{stage 2}).}
  \item{stop_s1}{Logical, indicating whether to stop after \ifelse{html}{\out{stage&nbsp;1}}{stage 1} 
  (due to \acronym{BE} or due to futility).}
  \item{stop_fut}{Logical, indicating whether study will be stopped after \ifelse{html}{\out{stage&nbsp;1}}{stage 1}
  due to futility.}
  \item{eRCI}{Exact repeated confidence interval for \ifelse{html}{\out{stage&nbsp;1}}{stage 1}.}
  \item{BE}{Logical, indicating whether \acronym{BE} can be concluded after 
  \ifelse{html}{\out{stage&nbsp;1}}{stage 1} or not.}
  \item{alpha_ssr}{Only applicable if study continues to \ifelse{html}{\out{stage&nbsp;2}}{stage 2}. Contains 
  alpha values for the two hypotheses required for sample size re-calculation.
  If \code{ssr.conditional = "no"} the result is equal to \code{alpha},
  otherwise it contains the conditional error rates for the standard combination
  test (in case of \code{max.comb.test = FALSE}) or maximum combination test (in
  case of \code{max.comb.test = TRUE}).}
  \item{GMR_ssr}{Only applicable if study continues to \ifelse{html}{\out{stage&nbsp;2}}{stage 2}. Contains the
  geometric mean ratio used for sample size re-calculation (accounts for
  adaptive planning step).}
  \item{targetpower_ssr}{Only applicable if study continues to \ifelse{html}{\out{stage&nbsp;2}}{stage 2}. 
  Contains the target power used for the sample size re-calculation (see also
  'Details').}
}
\author{
B. Lang
}
\references{
Patterson SD, Jones B. \emph{Bioequivalence and Statistics in Clinical Pharmacology.}\cr
Boca Raton: CRC Press; 2\ifelse{html}{\out{<sup>nd</sup>}}{\eqn{^{nd}}} edition 2017.

Maurer W, Jones B, Chen Y. \emph{Controlling the type 1 error rate in two-stage
sequential designs when testing for average bioequivalence.}\cr
Stat Med. 2018;1--21. \doi{10.1002/sim.7614}

Wassmer G, Brannath W. \emph{Group Sequential and Confirmatory Adaptive Designs
in Clinical Trials.}\cr
Springer 2016.
}
\seealso{
\code{\link{power.2stage.in, final.2stage.in}}
}
\examples{
# Example from Maurer et al
interim.2stage.in(GMR1 = exp(0.0424), CV1 = 0.3682, n1 = 20, max.n = 4000)
}
