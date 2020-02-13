\encoding{utf-8}
\name{interim.tsd.in}
\alias{interim.2stage.in}
\alias{interim.tsd.in}
\title{
Interim analysis of first stage data of 2-stage 2x2 crossover designs based on
the Inverse Normal method
}
\description{
  Following the design scheme according to \code{power.tsd.in} the function
  performs the interim analysis of the first stage data.
}
\usage{
interim.tsd.in(alpha, weight, max.comb.test = TRUE, targetpower = 0.8,
               GMR1, n1, CV1, df1 = NULL, SEM1 = NULL, theta1, theta2,
               GMR, usePE = FALSE, min.n2 = 4, max.n = Inf,
               fCpower = targetpower, fCrit = "CI", fClower, fCupper, fCNmax,
               ssr.conditional = c("error_power", "error", "no"),
               pmethod = c("exact", "nct", "shifted"))
}
\arguments{
  \item{alpha}{
    If one element is given, the overall one-sided significance level (not
    the adjusted level for \ifelse{html}{\out{stage&nbsp;1}}{stage 1}). In this
    case the adjusted alpha levels will be calcualted internally. If two
    elements are given, the argument refers to the two adjusted one-sided
    alpha levels for \ifelse{html}{\out{stage&nbsp;1}}{stage 1} and
    \ifelse{html}{\out{stage&nbsp;2}}{stage 2}, respectively.\cr
    If missing, defaults to \verb{0.05}.
  }
  \item{weight}{
    Pre-defined weight(s) of \ifelse{html}{\out{stage&nbsp;1}}{stage 1}, see
    'Details' for more information. Note that using the notation from
    Maurer et al, weight corresponds to information fraction, other literature
    may refer to sqrt(weight) as being the weight. \verb{weight} must either
    contain one element (in case of \verb{max.comb.test = FALSE}) or
    two elements (in case of \verb{max.comb.test = TRUE}).\cr
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
    Observed ratio of geometric means (T/R) of \ifelse{html}{\out{stage&nbsp;1}}{stage 1} data (use \emph{e.g.}, 0.95 for 95\%).
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
    \ifelse{html}{\out{stage&nbsp;1}}{stage 1} that can be specified in
    addition to \code{n1}.
  }
  \item{SEM1}{
    Optional; Standard error of the difference of means of
    \ifelse{html}{\out{stage&nbsp;1}}{stage 1} that can be specified in
    addition to \code{CV1}. Must be on additive scale (i.e. usually log-scale).
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
    If \verb{TRUE} the sample size re-estimation is done with the observed
    point estimate (\acronym{PE}) of the treatment difference in
    \ifelse{html}{\out{stage&nbsp;1}}{stage 1}.\cr
    Defaults to \verb{FALSE}.\cr
    Note: The power of stage 1 used for the futility inspection and calculation
    of the estimated conditional target power is always calculated with the
    planning value \code{GMR}.
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
    \ifelse{html}{\out{stage&nbsp;1}}{stage 1} is greater than or equal to
    \verb{fCpower}, then the study will be considered a failure.\cr\cr
    See \sQuote{Details} for more information on the choice of
    \verb{fCpower}.
  }
  \item{fCrit}{
    Futility criterion to use: \verb{"No"} (no futility criterion regarding
    observed point estimate, confidence interval and maximum sample size),
    \verb{"PE"} (observed point estimate of the geometric mean ratio from
    \ifelse{html}{\out{stage&nbsp;1}}{stage 1}),
    \verb{"CI"} (90\% confidence interval of the geometric mean ratio from
    \ifelse{html}{\out{stage&nbsp;1}}{stage 1}), \code{"Nmax"} (overall maximum
    sample size); or a combination thereof (concatenate abbreviations).
    Defaults to "CI".
  }
  \item{fClower}{
    Lower futility limit for the \acronym{PE} or \acronym{CI} of
    \ifelse{html}{\out{stage&nbsp;1}}{stage 1}.\cr
    If the \acronym{PE} or \acronym{CI} is completely outside of \verb{fClower}
    \ldots \verb{fCupper} the study is to be stopped due to futility (not \acronym{BE}).\cr
    May be missing. If \verb{"PE"} or \verb{"CI"} is specified within \verb{fCrit},
    the default will be set to 0.8 for \verb{fCrit = "PE"} or 0.95 for \verb{fCrit = "CI"}.
    If neither \verb{"PE"} nor \verb{"CI"} is specified within \verb{fCrit}, there
    will be no futility constraint regarding point estimate or confidence interval
    from \ifelse{html}{\out{stage&nbsp;1}}{stage 1} (regardless of any
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
    conditional error rates nor the estimated conditional target power for
    the second stage, \verb{"error"} uses conditional error rates for the
    second stage, and \verb{"error_power"} uses both conditional error rates
    and the estimated conditional target power for the second stage.\cr
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
    In contrast to \code{power.tsd.in} the default value here is \verb{"exact"}.
  }
}
\details{
  The observed values of stage 1 (e.g. \code{GMR1}, \code{n1}, \code{CV1}) may
  be obtained based on the first stage data via the usual ANOVA approach.\cr\cr
  The optional arguments \code{df1} and \code{SEM1} require a somewhat
  advanced knowledge (provided in the raw output from for example the software
  \acronym{SAS}, or may be obtained via \code{emmeans::emmeans}).
  However, it has the advantage that if there were missing data the exact
  degrees of freedom and standard error of the difference can be used,
  the former possibly being non-integer valued (e.g. if the
  Kenward-Roger method was used).\cr\cr
  The \code{weight} argument always refers to the first weight of a pair of
  weights. For example, in case of \code{max.comb.test = FALSE} the standard
  combination test requires two weights (w, 1-w) but only the first one, w,
  is required as input argument here because the second weight is
  automatically specified once the first is given. Similarly for
  \code{max.comb.test = TRUE}, w and w* need to be specified, which in turn
  define the two pairs of weights (w, 1-w) and (w*, 1-w*).\cr\cr
  If \verb{ssr.conditional = "error_power"}, the design scheme generally
  calculates the estimated conditional target power of the second stage and
  uses this value as desired target power in the sample size re-estimation process:\cr
  If \code{fCpower} > \code{targetpower}, then the conditional estimated
  target power may be negative. This does not seem sensible. Therefore, for such
  cases the desired target power for the sample size re-calculation will be set
  to \code{targetpower}, i.e. \code{ssr.conditional} will be set to \code{"error"}.\cr
  Also, if the futility criterion based on the power of stage 1 is met,
  then the conditional estimated target power will be negative. Thus, no further
  sample size calculation can be made. To acknowledge that this rule is nonbinding,
  for the purpose of calculating n2 the argument \code{ssr.conditional} is set
  to \code{"error"}.
}
\value{
  Returns an object of class \verb{"evaltsd"} with all the input arguments and results
  as components. As part of the input arguments a component \code{cval} is also
  presented, containing the critical values for
  \ifelse{html}{\out{stage&nbsp;1}}{stage 1} and 2 according to the
  input based on \code{alpha}, \code{weight} and \code{max.comb.test}.\cr
  The class \verb{"evaltsd"} has an S3 print method.\cr\cr
  The results are in the components:
  \item{p11}{Observed p-value for first hypothesis.}
  \item{p12}{Observed p-value for second hypothesis.}
  \item{z1}{z statistic value for first null hypothesis.}
  \item{z2}{z statistic value for second null hypothesis.}
  \item{RCI}{Repeated confidence interval for \ifelse{html}{\out{stage&nbsp;1}}{stage 1}. Corresponds to the usual CI with level alpha1.}
  \item{MEUE}{If the study stops, the median unbiased point estimate as estimate for the final adjusted
  geometric mean ratio after stage 1 (note that the value is identical to \code{GMR1}.)}
  \item{futility}{Three dimensional vector with either 0 or 1. The first
  component represents futility due to Power of first stage > \code{fCpower},
  the second futility due to \verb{CI} (or \verb{PE}) outside of
  \code{fClower ... fCupper}, the third futility due to
  \code{n1 + n2 > fCNmax}.\cr
  Note that the futility rules can be applied in a non-binding manner.}
  \item{CI90}{90\% Confidence interval for observed ratio of geometric means
  from \ifelse{html}{\out{stage&nbsp;1}}{stage 1}. If \verb{fCrit != "CI"} result will be \code{NULL}.}
  \item{Power Stage 1}{Calculated power of \ifelse{html}{\out{stage&nbsp;1}}{stage 1}.}
  \item{stop_s1}{Logical, indicating whether to stop after \ifelse{html}{\out{stage&nbsp;1}}{stage 1}
  (due to \acronym{BE} or due to futility).}
  \item{stop_fut}{Logical, indicating whether study is recommended to be stopped after
  \ifelse{html}{\out{stage&nbsp;1}}{stage 1}
  due to futility.}
  \item{stop_BE}{Logical, indicating whether \acronym{BE} could be concluded after
  \ifelse{html}{\out{stage&nbsp;1}}{stage 1} or not (regardless of any
  futility criterion).}
  \item{n2}{Required (total) sample size for \ifelse{html}{\out{stage&nbsp;2}}{stage 2}
  (will be zero if \acronym{BE} has been shown after
  \ifelse{html}{\out{stage&nbsp;1}}{stage 1}).}
  \item{alpha_ssr}{Only applicable if \acronym{BE} has not been shown after
  \ifelse{html}{\out{stage&nbsp;1}}{stage 1}. Contains
  alpha values for the two hypotheses required for sample size re-calculation.
  If \code{ssr.conditional = "no"} the result is equal to \code{alpha},
  otherwise it contains the conditional error rates for the standard combination
  test (in case of \code{max.comb.test = FALSE}) or maximum combination test (in
  case of \code{max.comb.test = TRUE}).}
  \item{GMR_ssr}{Only applicable if \acronym{BE} has not been shown after
  \ifelse{html}{\out{stage&nbsp;1}}{stage 1}. Contains the
  geometric mean ratio used for sample size re-calculation (accounts for
  adaptive planning step).}
  \item{targetpower_ssr}{Only applicable if \acronym{BE} has not been shown after
  \ifelse{html}{\out{stage&nbsp;1}}{stage 1}. Contains the target power used
  for the sample size re-calculation (see also 'Details').}
}
\author{
B. Lang
}
\references{
König F, Wolfsegger M, Jaki T, \enc{Schütz}{Schuetz} H, Wassmer G.\cr
\emph{Adaptive two-stage bioequivalence trials with early stopping and sample size re-estimation.}\cr
Vienna: 2014; 35\ifelse{html}{\out{<sup>th</sup>}}{\eqn{^{th}}} Annual Conference of the International Society for Clinical Biostatistics. Poster P1.2.88 \cr
\doi{10.13140/RG.2.1.5190.0967}.

Patterson SD, Jones B. \emph{Bioequivalence and Statistics in Clinical Pharmacology.}\cr
Boca Raton: CRC Press; 2\ifelse{html}{\out{<sup>nd</sup>}}{\eqn{^{nd}}} edition 2017.

Maurer W, Jones B, Chen Y. \emph{Controlling the type 1 error rate in two-stage
sequential designs when testing for average bioequivalence.}\cr
Stat Med. 2018; 37(10): 1587--1607. \doi{10.1002/sim.7614}.

Wassmer G, Brannath W. \emph{Group Sequential and Confirmatory Adaptive Designs
in Clinical Trials.}\cr
Springer 2016. \doi{10.1007/978-3-319-32562-0}.
}
\seealso{
\code{\link{power.tsd.in}, \link{final.tsd.in}}
}
\examples{
# Example from Maurer et al.
interim.tsd.in(GMR = 0.95, max.n = 4000,
               GMR1 = exp(0.0424), CV1 = 0.3682, n1 = 20)
# Example 2 from Potvin et al.
interim.tsd.in(GMR = 0.95, GMR1 = 1.0876, CV1 = 0.18213, n1 = 12,
               fCrit = "No", ssr.conditional = "no")
}
