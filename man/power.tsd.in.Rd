\encoding{utf-8}
\name{power.tsd.in}
\alias{power.2stage.in}
\alias{power.tsd.in}
\title{
Power calculation of adaptive 2-stage BE studies with 2x2 crossover design,
based on the Inverse Normal method
}
\description{
  The function calculates the \sQuote{empirical} power of 2-stage BE studies based
  on the Inverse-Normal combination method. The main design scheme is according
  to Maurer \emph{et al.} (Maximum Combination Test), but it may also be used 
  for other/modified designs, e.g. standard combination method, usage of the 
  observed treatment difference after \ifelse{html}{\out{stage&nbsp;1}}{stage 1} 
  in the sample size re-estimation step or different futility criteria.
}
\usage{
power.tsd.in(alpha, weight, max.comb.test = TRUE, n1, CV, targetpower = 0.8,
             theta0, theta1, theta2, GMR, usePE = FALSE, min.n2 = 4, max.n = Inf,
             fCpower = targetpower, fCrit = "CI", fClower, fCupper, fCNmax,
             ssr.conditional = c("error_power", "error", "no"),
             pmethod = c("nct", "exact", "shifted"), npct = c(0.05, 0.5, 0.95),
             nsims, setseed = TRUE, details = FALSE)
}
\arguments{
  \item{alpha}{
    If one element is given, the overall one-sided significance level. In this
    case the adjusted alpha levels will be calculated internally. If two
    elements are given, the argument refers to the two adjusted one-sided alpha 
    levels for \ifelse{html}{\out{stage&nbsp;1}}{stage 1} and
    \ifelse{html}{\out{stage&nbsp;2}}{stage 2}, respectively.\cr
    If missing, defaults to \verb{0.05}.
  }
  \item{weight}{
    Pre-defined weight(s) of \ifelse{html}{\out{stage&nbsp;1}}{stage 1}, see
    'Details' for more information.
    Note that using the notation from Maurer et al, weight corresponds to 
    information fraction, other literature may refer to sqrt(weight) as being
    the weight. \verb{weight} must either contain one element (in case of 
    \verb{max.comb.test = FALSE}) or two elements (in case of
    \verb{max.comb.test = TRUE}).\cr
    If missing, defaults to \verb{0.5} for \verb{max.comb.test = FALSE} and to
    \verb{c(0.5, 0.25)} for \verb{max.comb.test = TRUE}.
  }
  \item{max.comb.test}{
    Logical; if \verb{TRUE} (default) the maximum combination test will be
    used, otherwise the standard combination test.
  }
  \item{n1}{
    Sample size of \ifelse{html}{\out{stage&nbsp;1}}{stage 1}.
  }
  \item{CV}{
    Coefficient of variation of the intra-subject variability (use \emph{e.g.}, 0.3 for
    30\%).
  }
  \item{targetpower}{
    Desired (overall) target power to declare \acronym{BE} at the end of the trial.
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
    for \ifelse{html}{\out{stage&nbsp;1}}{stage 1} and sample size re-estimation
    for \ifelse{html}{\out{stage&nbsp;2}}{stage 2}. If missing, defaults to 0.95.
  }
  \item{usePE}{
    If \verb{TRUE} the sample size re-estimation is done with the observed 
    point estimate (\acronym{PE}) of the treatment difference in 
    \ifelse{html}{\out{stage&nbsp;1}}{stage 1}.\cr
    Defaults to \verb{FALSE}.\cr
    Note: The futility inspection via the power of stage 1 is always done with
    the planning value \code{GMR}.
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
    sample size); or a combination thereof (concatenate abbreviations; see
    \sQuote{Examples}). Defaults to "CI".
  }
  \item{fClower}{
    Lower futility limit for the \acronym{PE} or \acronym{CI} of
    \ifelse{html}{\out{stage&nbsp;1}}{stage 1}.\cr
    If the \acronym{PE} or \acronym{CI} is completely outside of \verb{fClower}
    \ldots \verb{fCupper} the study is stopped due to futility (not \acronym{BE}).\cr
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
    Defaults to \verb{"nct"} as it is a reasonable compromise between speed and
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
    the border of the acceptance range \verb{theta1} \ldots \verb{theta2}.
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
  But instead of simulating subject data, the statistics pe1, mse1 and pe2, mse2 
  are simulated via their associated distributions (Normal and
  \ifelse{html}{\out{&chi;<sup>2</sup>}}{\eqn{\chi^{2}}} distribution).\cr\cr
  The \code{weight} argument always refers to the first weight of a pair of
  weights. For example, in case of \code{max.comb.test = FALSE} the standard
  combination test requires two weights (w, 1-w) but only the first one, w,
  is required as input argument here because the second weight is
  automatically specified once the first is given. Similarly for 
  \code{max.comb.test = TRUE}, w and w* need to be specified, which in turn
  define the two pairs of weights (w, 1-w) and (w*, 1-w*).\cr\cr
  If \verb{ssr.conditional = "error_power"}, the design scheme generally 
  calculates the estimated conditional target power of the second stage and 
  uses this value as desired target power in the sample size re-estimation process. 
  If \code{fCpower} > \code{targetpower}, then the conditional power
  may actually be negative. This does not seem sensible. Therefore, for such
  cases the desired target power for the sample size re-calculation will be set 
  to \code{targetpower}.
}
\value{
  Returns an object of class \verb{"pwrtsd"} with all the input arguments and results
  as components. As part of the input arguments a component \code{cval} is also
  presented, containing the critical values for stage 1 and 2 according to the
  input based on \code{alpha}, \code{weight} and \code{max.comb.test}.\cr
  The class \verb{"pwrtsd"} has an S3 print method.\cr\cr
  The results are in the components:
  \item{pBE}{Fraction of studies found BE.}
  \item{pBE_s1}{Fraction of studies found BE in \ifelse{html}{\out{stage&nbsp;1}}{stage 1}.}
  \item{pct_stop_s1}{Percentage of studies stopped after
  \ifelse{html}{\out{stage&nbsp;1}}{stage 1} (due to BE or due to futility).}
  \item{pct_stop_fut}{Percentage of studies stopped after
  \ifelse{html}{\out{stage&nbsp;1}}{stage 1} due to futility.}
  \item{pct_s2}{Percentage of studies continuing to \ifelse{html}{\out{stage&nbsp;2}}{stage 2}.}
  \item{nmean}{Mean of n(total).}
  \item{nrange}{Range (min, max) of n(total).}
  \item{nperc}{Vector of percentiles of the distribution of n(total).}
}
\author{
B. Lang
}
\references{
König F, Wolfsegger M, Jaki T, Schütz H, Wassmer G. \cr
\emph{Adaptive two-stage bioequivalence trials with early stopping and sample size re-estimation.}\cr
Vienna: 2014; 35th Annual Conference of the International Society for Clinical Biostatistics. Poster P1.2.88 \cr 
\doi{10.13140/RG.2.1.5190.0967}.

Patterson SD, Jones B. \emph{Bioequivalence and Statistics in Clinical Pharmacology.}\cr
Boca Raton: CRC Press; 2\ifelse{html}{\out{<sup>nd</sup>}}{\eqn{^{nd}}} edition 2017.

Kieser M, Rauch G. \emph{Two-stage designs for cross-over bioequivalence trials.}\cr
Stat Med. 2015; 34(16):2403--16. \doi{10.1002/sim.6487}

Maurer W, Jones B, Chen Y. \emph{Controlling the type 1 error rate in two-stage
sequential designs when testing for average bioequivalence.}\cr
Stat Med. 2018;1--21. \doi{10.1002/sim.7614}

Wassmer G, Brannath W. \emph{Group Sequential and Confirmatory Adaptive Designs
in Clinical Trials.}\cr
Springer 2016. \doi{10.1007/978-3-319-32562-0}
}
\seealso{
\code{\link{interim.tsd.in}, \link{final.tsd.in}}
}
\examples{
# 12 subjects in stage 1, CV of 20\%, no futility rule, otherwise all defaults
# except too low number of sims due to CRAN policy for run-time of examples
# This setting resembles values in Table 6.8 and 6.10 in Patterson and Jones
# if nsims=1e5 (default)
power.tsd.in(n1 = 12, CV = 0.2, fCrit = "No", nsims=1e4)

# Note that actual replication would require the following additional
# arguments (extremely long run-time)
\dontrun{
power.tsd.in(n1 = 12, CV = 0.2, fCrit = "No", max.n = 4000, 
             pmethod = "exact", nsims = 1E6)}

# Table 8 in Maurer et al
power.tsd.in(n1 = 36, CV = 0.4, max.n = 4000)

# Same example as before but with additional futility criterion regarding
# maximum sample size (using the default 4*n1 as maximum bound)
power.tsd.in(n1 = 36, CV = 0.4, max.n = 4000, fCrit = c("CI", "Nmax"))
}
