\encoding{utf-8}
\name{power.2stage.KM}
\alias{power.2stage.KM}
\title{
Power calculation of adaptive 2-stage BE studies (2x2 crossover)
(variant of power.2stage to obtain the results of Karalis / Macheras)
}
\description{
This function calculates the empirical power of 2-stage BE studies according
to Potvin et.al. via simulations. The Potvin methods are modified as 
described in Karalis & Macheras to include a futility criterion Nmax 
and to do the power calculation steps and the sample size estimation step 
in the decision schemes with the mse (CV) \strong{and} the point estimate of T/R
from stage 1 evaluation.
}
\usage{
power.2stage.KM(method=c("C","B"), alpha0=0.05, alpha=c(0.0294,0.0294),
                n1, CV, targetpower=0.8, pmethod=c("nct","exact"),
                Nmax=150, theta0, theta1, theta2, npct=c(0.05, 0.5, 0.95), 
                nsims, setseed=TRUE, details=FALSE)
}
\arguments{
  \item{method}{
Decision schemes according to Potvin et.al.\cr
Default is "C" aka TSD in the paper of Karalis&Macheras if setting
\code{alpha=c(0.0294,0.0294)}.\cr
TSD-1 of the paper of Karalis can be obtained by choosing "C" but setting 
\code{alpha=c(0.028,0.028)}.\cr
TSD-2 of the paper of Karalis can be obtained by choosing "B" and setting  
\code{alpha=c(0.0294,0.0294)}.\cr
}
  \item{alpha0}{
Alpha value for the first step(s) in Potvin C aka TSD of Karalis&Macheras
or TSD-1 of the paper of Karalis, the power inspection and BE decision 
if power > targetpower.\cr
Defaults to 0.05.
}
  \item{alpha}{
Vector of the nominal alpha's for the two stages.\cr
Defaults to the Pocock alpha setting \code{alpha=c(0.0294,0.0294)}
as in TSD of the Karalis&Macheras paper.
}
  \item{n1}{
Sample size for stage 1.
}
  \item{CV}{
Coefficient of variation of the intra-subject variability as ratio.
}
  \item{targetpower}{
Power threshold in the first step of Potvin "C" and 
power to achieve in the sample size estimation step.
}
\item{pmethod}{
Power calculation method, also to be used in the sample size estimation for 
stage 2.\cr 
Implemented are "nct" = approximate calculations via non-central t-distribution 
and "exact" = exact calculations via Owen's Q.\cr
Defaults to "nct" for speed reasons in the sample size estimation step.
}
  \item{Nmax}{
Futility criterion. If set to a finite value all studies simulated in which a 
sample size >Nmax is obtained will be regarded as BE=FAIL.\cr
Defaults to 150, recommended in the Karalis & Macheras paper.\cr
Set this argument to \code{Inf}, to work without that futility criterion.
}
  \item{theta0}{
True ratio of T/R for simulating. Defaults to 0.95.
}
  \item{theta1}{
Lower bioequivalence limit. Defaults to 0.8.
}
  \item{theta2}{
Upper bioequivalence limit. Defaults to 1.25.
}
  \item{npct}{
Percentiles to be used for the presentation of the distribution of n(total)=n1+n2.\cr
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
The calculations follow in principle the simulations as described in Potvin 
et. al.\cr
The underlying subject data are assumed to be evaluated after log-transformation.
But instead of simulating subject data the statistics pe1, mse1 and pe2, SS2 
are simulated via their associated distributions 
(normal and chi-squared distribution's).\cr\cr
In contrast to Potvin et.al. the power calculation steps as well as the sample 
size adaption step of the decision schemes are done using the the mse (CV) 
\bold{and} the point estimate from stage 1.\cr
This resembles the methods described in Karalis&Macheras and Karalis.  
}
\value{
Returns an object of class "pwrtsd" with all the input arguments and results 
as components.\cr
The class "pwrtsd" has a S3 print method.\cr
The results are in the components:
\item{pBE}{Contains the ratio of studies found BE.}
\item{pBE_s1}{Ratio of studies found BE in stage 1.}
\item{pct_s2}{Percentage of studies continuing to stage 2.}
\item{nmean}{Mean of n(total).}
\item{nrange}{Range (min, max) of n(total).}
\item{nperc}{Percentiles of the distribution of n(total).}
\item{ntable}{Object of class "table" summarizing the discrete distribution of 
n(total) via its distinct values and counts of occurences of these values.\cr
This component is only given back if \code{is.finite(Nmax)}.}
}
\references{
Potvin D et al.\cr
\emph{Sequential design approaches for bioequivalence studies with 
crossover designs}\cr
Pharm Stat 7(4), 245-62 (2008) \doi{10.1002/pst.294}

Karalis V and P Macheras\cr
\emph{An Insight into the Properties of a Two-Stage Design in Bioequivalence Studies}\cr
Pharm. Res. 30(7) 1824-35 (2013), Epub 2013 April 9. \doi{10.1007/s11095-013-1026-3}

Karalis V\cr
\emph{The role of the upper sample size limit in two-stage bioequivalence designs}\cr
Int J Pharm. 2013 Nov 1;456(1):87-94. Epub 2013 Aug 14. \doi{10.1016/j.ijpharm.2013.08.013}

Fuglsang A\cr
\emph{Futility Rules in Bioequivalence Trials with Sequential Designs}\cr
AAPS J. 2014 Jan;16(1):79-82. Epub 2013 Nov 12. \doi{10.1208/s12248-013-9540-0}
}
\author{
D. Labes
}
\examples{
# using all the defaults
# but too low number of sims to complain with the CRAN policy: 
# "check time only a few seconds per example"
# minimum number of sims should be 1E5 for power, 1E6 sims for 'alpha'
power.2stage.KM(n1=16, CV=0.2, nsims=1E4)
# 3.3 secs if nsims=1E5
}
