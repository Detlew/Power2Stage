\encoding{utf-8}
\name{power.2stage.in}
\alias{power.2stage.in}
\title{
Power calculation of adaptive 2-stage BE studies with 2x2 crossover design,
based on the Inverse-Normal method
}
\description{
The function calculates the empirical power of 2-stage BE studies based
on the Inverse-Normal combination method. The main design scheme is according
to Maurer et al, but it may also be used for modified designs to allow for
usage of the observed treatment difference after stage 1 in the sample size
re-estimation step, or a different power threshold to decide on futility after
stage 1.
}
\usage{
power.2stage.in(alpha, weight, max.comb.test = TRUE, n1, CV, targetpower = 0.8,
                power.threshold = targetpower, theta0, theta1, theta2, GMR, 
                usePE = FALSE, min.n2 = 4, max.n = Inf, ssr.conditional = TRUE, 
                fCrit = NULL, fClower, fCupper, fCNmax, 
                pmethod = c("nct", "exact", "shifted"), 
                npct = c(0.05, 0.5, 0.95), 
                nsims, setseed = TRUE, details = FALSE)
}
\arguments{
  \item{alpha}{
If of length 1, overall one-sided significance level. If of length 2, adjusted
one-sided alpha levels for stage 1 and 2.\cr
If missing, the default value is set to \code{alpha = 0.05}. 
}
  \item{weight}{
Weight(s) of stage 1 (information fraction). This argument will be ignored if
\code{alpha} is of length 2. Otherwise, \code{weight} must either be of
length 1 (in case of \code{max.comb.test = FALSE}) or of length 2 (in case of
\code{max.comb.test = TRUE}).\cr
If missing, the default is 0.5 for \code{max.comb.test = FALSE} or 
\code{c(0.5, 0.25)} for \code{max.comb.test = TRUE}.
}
  \item{max.comb.test}{
Logical; if \code{TRUE} the maximum combination test will be used, otherwise
the standard combination test.
}
  \item{n1}{
Sample size for stage 1.
}
  \item{CV}{
Coefficient of variation of the intra-subject variability (use e.g. 0.3 for
30\%).
}
  \item{targetpower}{
Desired target power to declare BE at the end of the trial.
}
  \item{power.threshold}{
Threshold for power monitoring step to decide on futility for cases where
BE has not been achieved after stage 1: If BE has not been achieved after
stage 1 and power for stage 1 is greater than \code{power.threshold}, then the 
study will be considered a failure.\cr
NB: Set value to 1 in order to deactivate this futility rule.
}
  \item{theta0}{
Assumed ratio of geometric means (T/R) for simulations.\cr
If missing, the default value is equal to argument \code{GMR}.
}
  \item{theta1}{
Lower bioequivalence limit. Defaults to 0.8.
}
  \item{theta2}{
Upper bioequivalence limit. Defaults to 1.25.
}
  \item{GMR}{
Assumed ratio of geometric means (T/R) to be used in power calculation for 
stage 1 and sample size re-estimation step for stage 2.
}
  \item{usePE}{
If \code{TRUE} the sample size re-estimation step is done with the observed 
point estimate (PE) for the treatment difference from stage 1.\cr
Defaults to \code{FALSE}.
}
  \item{min.n2}{
Minimum sample size of stage 2. Defaults to 4.\cr
If the sample size re-estimation step gives a sample size for stage 2 less
than \code{min.n2}, then \code{min.n2} will be used for stage 2.
}
  \item{max.n}{
Maximum overall sample size (stage 1 + stage 2).\cr
This is \emph{not} a futility criterion regarding the maximum sample size! If 
max.n is set to a finite value and the sample size re-estimation step gives a 
sample size for stage 2, n2, such that n1 + n2 > max.n, then the sample size 
for stage 2 will be set to n2 = max.n - n1.\cr
Defaults to Inf, i.e. no constraint on the re-estimated sample size.
}
  \item{ssr.conditional}{
Logical; if \code{TRUE}, the sample size re-estimation step uses conditional
error rates and conditional power.\cr
Defaults to \code{TRUE}.
}
  \item{fCrit}{
Futility criterion to use: \code{"PE"}, \code{"CI"} or \code{"Nmax"}; 
or a combination thereof (concatenate abbreviations).
}
  \item{fClower}{
Lower futility limit for the PE or CI of stage 1.\cr
If the PE or CI is completely outside of \code{fClower} ... \code{fCupper} 
the study is stopped due to futility (no BE).\cr
May be missing. If "PE" or "CI" is specified within \code{fCrit}, the default 
will be set to 0.8 for \code{fCrit = "PE"} or 0.95 for \code{fCrit = "CI"}.
If neither "PE" nor "CI" is specified within \code{fCrit}, there will be
no futility constraint regarding point estimate or confidence interval from
stage 1 (regardless of any input for \code{fCupper}).
}
  \item{fCupper}{
Upper futility limit for the PE or CI of stage 1.\cr
Analogous to \code{fClower}: will be set to \code{1/fClower} if missing.
}
  \item{fCNmax}{
Futility criterion regarding maximum sample size. If the determined sample size
for stage 2, n2, is such that n1 + n2 > fCNmax, the study will not continue to
stage 2 and will be stopped due to futility (no BE).\cr
If Nmax is specified within \code{fCrit} and argument \code{fCNmax} is missing,
the value will be set to fCNmax = 4*n1. If Nmax is not specified within
\code{fCrit}, then there will be no futility constraint regarding maximum
sample size (regardless of any input for \code{fCNmax}).
}
  \item{pmethod}{
Power calculation method, also to be used in the sample size estimation for 
stage 2.\cr 
Implemented are \code{"nct"} = approximate calculations via non-central 
t-distribution and \code{"exact"} = exact calculations via Owen's Q as well as
the approximation via shifted central t-distribution, as used by Potvin et al.\cr
Defaults to "nct" as reasonable compromise between speed and accuracy in the 
sample size estimation step.
}
  \item{npct}{
Percentiles to be used for the presentation of the distribution of 
n(total)=n1+n2.\cr
Defaults to \code{c(0.05, 0.5, 0.95)} to obtain the 5\% and 95\% percentiles
and the median.
}
  \item{nsims}{
Number of studies to simulate.\cr 
If missing is set to 1E+05 = 100 000 or to 1E+06 = 1 mio if you are calculating 
type 1 error, i.e. with \code{theta0} at border or outside acceptance range 
\code{theta1} ... \code{theta2}.
}
  \item{setseed}{
Simulations are dependent on the starting point of the (pseudo) random number 
generator. To avoid differences in power for different runs a 
\code{set.seed(1234567)} is issued if \code{setseed=TRUE}, the default.\cr
Set this argument to \code{FALSE} to view the variation in power between 
different runs.
}
  \item{details}{
If set to \code{TRUE} the function prints the results of time measurements
of the simulation steps. Default is \code{FALSE}.
}
}
\details{
The underlying subject data are assumed to be evaluated after log-transformation.
But instead of simulating subject data, the statistics pe1, mse1 and pe2, mse2 
are simulated via their associated distributions (normal and chi-squared 
distributions).
}
\value{
Returns an object of class "pwrtsd" with all the input arguments and results 
as components.\cr
The class "pwrtsd" has an S3 print method.\cr
The results are in the components:
\item{pBE}{Contains the ratio of studies found BE.}
\item{pBE_s1}{Ratio of studies found BE in stage 1.}
\item{pct_stop_s1}{Percentage of studies stopped after stage 1 (due to BE or
due to futility).}
\item{pct_stop_fut}{Percentage of studies stopped after stage 1 due to
futility.}
\item{pct_s2}{Percentage of studies continuing to stage 2.}
\item{nmean}{Mean of n(total).}
\item{nrange}{Range (min, max) of n(total).}
\item{nperc}{Vector of percentiles of the distribution of n(total).}
}
\references{
Patterson S, Jones B\cr
\emph{Bioequivalence and Statistics in Clinical Pharmacology}\cr
Chapman & Hall/CRC, Boca Raton, Second Edition 2017.

Maurer W, Jones B, Chen Y\cr
(To be Submitted) \emph{Controlling the type 1 error rate in two-stage 
sequential designs when testing for average bioequivalence}\cr
Statistics in Medicine.

Kieser M, Rauch G\cr
\emph{Two-stage designs for cross-over bioequivalence trials}\cr
Stat Med. 2015 Jul 20;34(16):2403-16\cr
Epub 2015 Mar 24 \doi{10.1002/sim.6487}
}
\author{
B. Lang
}
\examples{
# Using all defaults, 12 subjects in stage 1 and a CV of 20\%
power.2stage.in(n1 = 12, CV = 0.2)

# This setting resembles values in Table 6.8 and 6.10 in Patterson and Jones. 
# Note that actual replication would require the following additional 
# arguments (long run-time)
\dontrun{
power.2stage.in(n1 = 12, CV = 0.2, max.n = 4000, pmethod = "exact", nsims = 1E6)}

# Example as above, but using futility constraints (the default ones) 
# regarding PE and Nmax:
power.2stage.in(n1 = 12, CV = 0.2, fCrit = c("PE", "Nmax"))
}