\encoding{utf-8}
\name{power.2stage.in}
\alias{power.2stage.in}
\title{
Power calculation of adaptive 2-stage BE studies with 2x2 crossover design,
based on the Inverse-Normal method
}
\description{
  The function calculates the \sQuote{empiric} power of 2-stage BE studies based
  on the Inverse-Normal combination method. The main design scheme is according
  to Maurer \emph{et al.}, but it may also be used for modified designs to
  allow for usage of the observed treatment difference after
  \ifelse{html}{\out{stage&nbsp;1}}{stage~1} in the sample size re-estimation
  step, or a different power threshold to decide on futility after
  \ifelse{html}{\out{stage&nbsp;1}}{stage~1}.
}
\usage{
power.2stage.in(alpha, weight, max.comb.test = TRUE, n1, CV, targetpower = 0.8,
                power.threshold = targetpower, theta0, theta1, theta2, GMR,
                usePE = FALSE, min.n2 = 4, max.n = Inf, ssr.conditional = TRUE,
                fCrit = NULL, fClower, fCupper, fCNmax,
                pmethod = c("nct", "exact", "shifted"), npct = c(0.05, 0.5, 0.95),
                nsims, setseed = TRUE, details = FALSE)
}
\arguments{
  \item{alpha}{
    If one element is given, the overall one-sided significance level. If two
    element are given, the adjusted one-sided alpha levels for
    \ifelse{html}{\out{stage&nbsp;1}}{stage~1} and
    \ifelse{html}{\out{stage&nbsp;2}}{stage~21}, respectively.\cr
    If missing, defaults to \verb{0.05}.
  }
  \item{weight}{
    Weight(s) of \ifelse{html}{\out{stage&nbsp;1}}{stage~1} (information fraction).
    This argument will be ignored if \verb{alpha} contains two elements.
    Otherwise, \verb{weight} must either be either contain one element (in case
    of \verb{max.comb.test = FALSE}) or two elements (in case of
    \verb{max.comb.test = TRUE}).\cr
    If missing, defaults to \verb{0.5} (for \verb{max.comb.test = FALSE}) or
    \verb{c(0.5, 0.25)} (for \verb{max.comb.test = TRUE}).
  }
  \item{max.comb.test}{
    Logical; if \verb{TRUE} (default) the maximum combination test will be
    used, otherwise the standard combination test.
  }
  \item{n1}{
    Sample size of \ifelse{html}{\out{stage&nbsp;1}}{stage~1}.
  }
  \item{CV}{
    Coefficient of variation of the intra-subject variability (use \emph{e.g.}, 0.3 for
    30\%).
  }
  \item{targetpower}{
    Desired target power to declare \acronym{BE} at the end of the trial.
  }
  \item{power.threshold}{
    Threshold for power monitoring step to decide on futility for cases where
    \acronym{BE} has not been achieved after
    \ifelse{html}{\out{stage&nbsp;1}}{stage~1}: If \acronym{BE} has not been
    achieved after \ifelse{html}{\out{stage&nbsp;1}}{stage~1} and power for
    \ifelse{html}{\out{stage&nbsp;1}}{stage~1} is greater than or equal to
    \verb{power.threshold}, then the study will be considered a failure.\cr\cr
    See \sQuote{Details} for more information on the choice of
    \verb{power.threshold}.
  }
  \item{theta0}{
    Assumed ratio of geometric means (T/R) for simulations. If missing,
    defaults to \verb{GMR}.
  }
  \item{theta1}{
    Lower bioequivalence limit. Defaults to 0.8.
  }
  \item{theta2}{
    Upper bioequivalence limit. Defaults to 1.25.
  }
  \item{GMR}{
    Assumed ratio of geometric means (T/R) to be used in power calculation
    for \ifelse{html}{\out{stage&nbsp;1}}{stage~1} and sample size re-estimation
    for \ifelse{html}{\out{stage&nbsp;2}}{stage~2}.
  }
  \item{usePE}{
    If \verb{TRUE} the sample size re-estimation step is done with the
    observed point estimate (\acronym{PE}) of the treatment difference in
    \ifelse{html}{\out{stage&nbsp;1}}{stage~1}.\cr
    Defaults to \verb{FALSE}.
  }
  \item{min.n2}{
    Minimum sample size of \ifelse{html}{\out{stage&nbsp;2}}{stage~2}. Defaults to 4.\cr
    If the sample size re-estimation step gives a sample size for
    \ifelse{html}{\out{stage&nbsp;2}}{stage~2} less than \verb{min.n2}, then
    \verb{min.n2} will be used for \ifelse{html}{\out{stage&nbsp;2}}{stage~2}.
  }
  \item{max.n}{
    Maximum overall sample size \ifelse{html}{\out{stage&nbsp;1}}{stage~1} +
    \ifelse{html}{\out{stage&nbsp;2}}{stage~2}).\cr
    This is \emph{not} a futility criterion regarding the maximum sample size! If
    \verb{max.n} is set to a finite value and the sample size re-estimation gives a
    sample size for \ifelse{html}{\out{stage&nbsp;2}}{stage~2} (\verb{n2}), such
    that \verb{n1 + n2 > max.n}, then the sample size for
    \ifelse{html}{\out{stage&nbsp;2}}{stage~2} will be set to \verb{n2 = max.n - n1}.\cr
    Defaults to \verb{Inf}, \emph{i.e.}, no constraint on the re-estimated sample size.
  }
  \item{ssr.conditional}{
    Logical; if \verb{TRUE}, the sample size re-estimation step uses conditional
    error rates and conditional power. Defaults to \verb{TRUE}.\cr
    See \sQuote{Details}.
  }
  \item{fCrit}{
    Futility criterion to use: \verb{"PE"} (observed point estimate of the
    geometric mean ratio from \ifelse{html}{\out{stage&nbsp;1}}{stage~1}),
    \verb{"CI"} (90\% confidence interval of the geometric mean ratio from
    \ifelse{html}{\out{stage&nbsp;1}}{stage~1}) or \code{"Nmax"} (overall maximum
    sample size); or a combination thereof (concatenate abbreviations; see
    \sQuote{Examples}).
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
    from \ifelse{html}{\out{stage&nbsp;1}}{stage~1} (regardless of any specification
    of \verb{fCupper}).
  }
  \item{fCupper}{
    Upper futility limit for the \acronym{PE} or \acronym{CI} of
    \ifelse{html}{\out{stage&nbsp;1}}{stage~1}.\cr
    Analogous to \verb{fClower}: Will be set to \verb{1/fClower} if missing.
  }
  \item{fCNmax}{
    Futility criterion regarding maximum sample size. If the determined sample size
    for \ifelse{html}{\out{stage&nbsp;2}}{stage~2} (\verb{n2}), is such that
    \verb{n1 + n2 > fCNmax}, the study will not continue to
    \ifelse{html}{\out{stage&nbsp;2}}{stage~2} and stopped due to futility (not
    \acronym{BE}).\cr
    If \verb{Nmax} is specified within \verb{fCrit} and argument \verb{fCNmax}
    is missing, the value will be set to \verb{fCNmax = 4*n1}. If \verb{Nmax} is
    not specified within \verb{fCrit}, then there will be no futility constraint
    regarding maximum sample size (regardless of any specification of \verb{fCNmax}).
  }
  \item{pmethod}{
    Power calculation method, also to be used in the sample size estimation for
    \ifelse{html}{\out{stage&nbsp;2}}{stage~2}.\cr
    Implemented are \verb{"nct"} (approximate calculations via non-central
    \emph{t}-distribution, \verb{"exact"} (exact calculations via Owen\enc{’}{'}s Q),
    and \verb{"shifted"} (approximate calculation via shifted central \emph{t}-distribution
    like in the paper of Potvin \emph{et al.}\cr
    Defaults to \verb{"nct"} as a reasonable compromise between speed and
    accuracy in the sample size estimation step.
  }
  \item{npct}{
    Percentiles to be used for the presentation of the distribution of
    \verb{n(total)=n1+n2}.\cr
    Defaults to \verb{c(0.05, 0.5, 0.95)} to obtain the 5\% and 95\% percentiles
    and the median.
  }
  \item{nsims}{
    Number of studies to simulate.\cr
    If missing, \verb{nsims} is set to 1E+05 = 100,000 or to 1E+06 = 1 Mio if
    estimating the empiric Type I Error (\verb{'alpha'}), \emph{i.e.}, with \verb{theta0} at
    the border or outside the acceptance range \verb{theta1} \ldots \verb{theta2}.
  }
  \item{setseed}{
    Simulations are dependent on the starting point of the (pseudo) random number
    generator. To avoid differences in power for different runs a
    \code{set.seed(1234567)} is issued if \verb{setseed=TRUE}, the default.\cr
    Set this argument to \verb{FALSE} to view the variation in power between
    different runs.
  }
  \item{details}{
  If set to \verb{TRUE} the function prints the results of time measurements
  of the simulation steps. Default to \verb{FALSE}.
  }
}
\details{
  The underlying subject data are assumed to be evaluated after log-transformation.
  But instead of simulating subject data, the statistics pe1, mse1 and pe2, SS2 are
  simulated via their associated distributions (normal and
  \ifelse{html}{\out{&chi;<sup>2</sup>}}{\eqn{\chi^{2}}} distributions).\cr\cr
  If \verb{ssr.conditional = TRUE}, the design scheme generally calculates the
  conditional power of the second stage and uses this value as desired target
  power in the sample size re-estimation process. There are two links between
  the conditional power and the argument \verb{power.threshold}:
  \itemize{
  \item If \code{power.threshold} > \code{targetpower}, then the conditional power
  may actually be negative. This seems not to be sensible. Therefore, for such
  cases the desired target power for sample size re-calculation will be set to
  \code{targetpower}.
  \item The futility criterion \eqn{power for stage 1 \ge } \code{power.threshold}
  for cases with not yet BE after stage 1 may equivalently be expressed in terms of
  the conditional power of the second stage: The criterion
  \deqn{conditional power \le x} is equivalent to
  \deqn{power for stage 1 \ge 1 - \beta/(1 - x),} where
  \eqn{\beta = 1 - }\code{targetpower} and \eqn{x in [0, 1)}.
  Therefore, \code{power.threshold} may be chosen such that futility is based
  on the conditional power being less than some value \eqn{x}.
  }
}
\value{
  Returns an object of class \verb{"pwrtsd"} with all the input arguments and results
  as components. As part of the input arguments a component \code{cval} is also
  presented, containing the critical values for stage 1 and 2 according to the
  input based on \code{alpha}, \code{weight} and \code{max.comb.test}.\cr
  The class \verb{"pwrtsd"} has an S3 print method.\cr\cr
  The results are in the components:
  \item{pBE}{Fraction of studies found BE.}
  \item{pBE_s1}{Fraction of studies found BE in \ifelse{html}{\out{stage&nbsp;1}}{stage~1}.}
  \item{pct_stop_s1}{Percentage of studies stopped after
  \ifelse{html}{\out{stage&nbsp;1}}{stage~1} (due to BE ordue to futility).}
  \item{pct_stop_fut}{Percentage of studies stopped after
  \ifelse{html}{\out{stage&nbsp;1}}{stage~1} due to futility.}
  \item{pct_s2}{Percentage of studies continuing to \ifelse{html}{\out{stage&nbsp;2}}{stage~2}.}
  \item{nmean}{Mean of n(total).}
  \item{nrange}{Range (min, max) of n(total).}
  \item{nperc}{Vector of percentiles of the distribution of n(total).}
}
\author{
B. Lang
}
\references{
Patterson SD, Jones B. \emph{Bioequivalence and Statistics in Clinical Pharmacology.}\cr
Boca Raton: CRC Press; 2\ifelse{html}{\out{<sup>nd</sup>}}{\eqn{^{nd}}} edition 2016.

Maurer W, Jones B, Chen Y. \emph{Controlling the type 1 error rate in two-stage
sequential designs when testing for average bioequivalence.}\cr
Stat Med. (to be submitted).

Kieser M, Rauch G. \emph{Two-stage designs for cross-over bioequivalence trials.}\cr
Stat Med. 2015; 34(16):2403--16. \doi{10.1002/sim.6487}
}

\examples{
# Using all defaults, 12 subjects in stage 1 and a CV of 20\%
power.2stage.in(n1 = 12, CV = 0.2)

# This setting resembles values in Table 6.8 and 6.10 in Patterson and Jones.
# Note that actual replication would require the following additional
# arguments (extremely long run-time)
\dontrun{
power.2stage.in(n1 = 12, CV = 0.2, max.n = 4000, pmethod = "exact", nsims = 1E6)}
#
# Example as above, but using futility constraints (the default ones)
# regarding PE and Nmax:
power.2stage.in(n1 = 12, CV = 0.2, fCrit = c("PE", "Nmax"))
}
