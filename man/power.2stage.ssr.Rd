\encoding{utf-8}
\name{power.2stage.ssr}
\alias{power.2stage.ssr}
\title{
Power of 2-stage BE studies with 2x2 crossover design with interim
sample size re-estimation
}
\description{
This function calculates the empirical power (via simulations) of 2-stage BE 
studies with interim sample size re-estimation (no BE decision after stage 1).
The sample size re-estimation may be done blinded or not blinded.
}
\usage{
power.2stage.ssr(alpha = 0.05, n1, GMR, CV, targetpower = 0.8, 
                pmethod = c("nct", "exact", "shifted", "ls"), blind = FALSE, 
                usePE = FALSE, min.n = 0, max.n = Inf, theta0, theta1, theta2, 
                npct = c(0.05, 0.5, 0.95), nsims, setseed = TRUE, details = FALSE)
}
\arguments{
  \item{alpha}{
Nominal type I error. Has to be adapted in case of alpha-inflation.
}
  \item{n1}{
Sample size of stage 1.
}
  \item{GMR}{
Ratio T/R to be used in the sample size re-estimation.\cr
Defaults to 0.95 if missing.
}
  \item{CV}{
Coefficient of variation of the intra-subject variability as ratio.\cr
Anticipated population value.
}
  \item{targetpower}{
Power to achieve in the sample size re-estimation step.
}
  \item{pmethod}{
Power calculation method to be used in the sample size re-estimation for stage 2.\cr 
Implemented are \code{"nct"} = approximate calculations via non-central 
t-distribution and \code{"exact"} = exact calculations via Owen's Q as well as
the approximation via shifted central t-distribution.\cr
Also implemented is the large sample approximation as used in the references.\cr
Defaults to \code{"nct"}.
}
  \item{blind}{
If \code{TRUE} the blinded estimate of the intra-subject variance, i.e. the 
estimate from the period differences, is used in sample size estimation.\cr
If \code{FALSE} the usual mse from stage 1 is used.\cr
Defaults to \code{FALSE} since most BE studies are done as open.
}
  \item{usePE}{
If \code{TRUE} the point estimate from the interim analysis is used in the 
sample site re-estimation. Defaults to \code{FALSE}.\cr\cr
\code{usePE = TRUE} doesn't make sense if \code{blind = TRUE}. In that case the
function issues a warning and \code{usePE} is reset to \code{usePE = FALSE}.
}
  \item{min.n}{
If \code{min.n>n1}, the re-estimated sample size (N) is set to max(min.n,N).\cr
If \code{min.n=0}, the default, no minimal sample size is applied.
}
  \item{max.n}{
If max.n is set to a finite value the re-estimated sample size (N) is set to 
min(max.n,N).\cr
Defaults to Inf which is equivalent to not constrain the re-estimated sample size.\cr
Attention! max.n is here \bold{not} a futility criterion like Nmax in other 
functions of this package.
}
  \item{theta0}{
True ratio of T/R for simulating. Defaults to the \code{GMR} argument if missing.
}
  \item{theta1}{
Lower bioequivalence limit. Defaults to 0.8.
}
  \item{theta2}{
Upper bioequivalence limit. Defaults to 1.25.
}
  \item{npct}{
p-values to be used for the percentiles of the distribution of n(total)=n1+n2.\cr
Defaults to \code{c(0.05, 0.5, 0.95)} to obtain the 5\% and 95\% percentiles
and the median.
}
  \item{nsims}{
Number of studies to simulate.\cr 
If missing \code{nsims} is set to 1E+05 = 100 000 or to 1E+06 = 1 Mio if you are 
calculating 'alpha', i.e. with \code{theta0} at border or outside acceptance range 
\code{theta1} ... \code{theta2}.
}
  \item{setseed}{
Simulations are dependent on the starting point of the (pseudo) random number 
generator. To avoid differences in power for different runs a \code{set.seed(1234567)} 
is issued if \code{setseed=TRUE}, the default.\cr
Set this argument to \code{FALSE} to view the variation in power between 
different runs.
}
  \item{details}{
If set to \code{TRUE} the function prints the results of time measurements
of the simulation steps. Default is \code{FALSE}.
}
}
\details{
The calculations follow in principle the simulations as described in 
Potvin et al. but with no power monitoring step(s) and no BE decision in 
stage 1 evaluation.\cr
The underlying subject data are assumed to be evaluated after log-transformation.
But instead of simulating subject data the statistics pe1, mse1 and pe2, SS2 are
simulated via their associated distributions (normal and chi-squared 
distri's).
}
\value{
Returns an object of class "pwrtsd" with all the input arguments and results 
as components.\cr
The class "pwrtsd" has a S3 print method.\cr
The results are in the components:
\item{pBE}{Contains the ratio of studies found BE.}
\item{pct_s2}{Percentage of studies continuing to stage 2.}
\item{nmean}{Mean of n(total).}
\item{nrange}{Range (min, max) of n(total).}
\item{nperc}{Percentiles of the distribution of n(total).}
\item{ntable}{Object of class "table" summarizing the discrete distribution of 
n(total) via its unique values and counts of occurences of these values.\cr
ntable is only given back if \code{usePE = FALSE}}
}
\references{
Golkowski D, Friede T, Kieser M\cr
\emph{Blinded sample size re-estimation in crossover bioequivalence trials}\cr
Pharm Stat. 2014 May-Jun;13(3):157-62. Epub 2014 Apr 9. \doi{10.1002/pst.1617}

Jones B, Kenward MG\cr
\emph{Design and Analysis of Cross-Over Trials}\cr
Third edition, Chapman and Hall/CRC, Boca Raton, Oct 2014\cr
Chapter 12

Potvin D et al.\cr
\emph{Sequential design approaches for bioequivalence studies with 
crossover designs}\cr
Pharm Stat 7(4), 245-62 (2008) \doi{10.1002/pst.294}
}
\author{
D. Labes
}
\note{
The computation time is in the magnitude of a few seconds for 100 000 sims on my 
machine (Intel core i7 2.5 GHz, 12GB RAM) if the nct approximation is used. 
Thus be a bit patient if you simulate for 'alpha' with 1 Mio sims.\cr
Using the crude \code{pmethod="ls"} on the other hand results in a nearly immediate 
sample size re-estimation.
}
\seealso{
\code{\link{power.2stage}} for 2-stage studies with interim BE decision.
}
\examples{
# much too low number of sims to complain with CRAN policy about examples run-time 
# minimum number of sims should be 1E5 for 'power', 1E6 sims for 'alpha'
power.2stage.ssr(alpha=0.05, n1=10, CV=0.245, GMR=1, pmethod="ls", blind=TRUE,
                 theta0=1.25, nsims=1E4)
# computation time ~ 5 secs if 1E6 sims were used with pmethod="ls"
# should give the maximum alpha-inflation ~ 0.0723 observed in Golkowski et al.
# for blinded sample size re-estimation if nsims=1E6 is used.
}
