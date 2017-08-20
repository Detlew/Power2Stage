\encoding{utf-8}
\name{power.2stage}
\alias{power.2stage}
\title{
Power calculation of adaptive 2-stage BE studies with 2x2 crossover design
}
\description{
This function calculates the empirical power of 2-stage BE studies according
to Potvin et al. via simulations. The Potvin methods are modified to include
a futility criterion Nmax and to allow the sample size estimation step to be 
done with point estimate & mse of stage 1.
}
\usage{
power.2stage(method = c("B", "C", "B0"), alpha0 = 0.05, alpha = c(0.0294, 0.0294), 
             n1, GMR, CV, targetpower = 0.8, pmethod = c("nct", "exact", "shifted"),
             usePE = FALSE, Nmax = Inf, min.n2=0, theta0, theta1, theta2,  
             npct = c(0.05, 0.5, 0.95), nsims, setseed = TRUE, details = FALSE)
}
\arguments{
  \item{method}{
Decision schemes according to Potvin et.al. Default is "B". \cr
Potvin method D can be obtained by choosing "C" but setting 
\code{alpha=c(0.028,0.028)}\cr.
\code{method="B0"} uses the decision scheme of so-called "MSDBE" which differs
from B/E in case of different alpha's w.r.t. power monitoring and BE decision
in case of power >= target power.
}
  \item{alpha0}{
Alpha value for the first step(s) in Potvin "C"", the power inspection
and BE decision if power > targetpower.\cr
Defaults to 0.05.
}
  \item{alpha}{
Vector of the nominal alpha's for the two stages.\cr
Defaults to the Pocock alpha setting \code{alpha=c(0.0294,0.0294)}.
}
  \item{n1}{
Sample size for stage 1.
}
  \item{GMR}{
Ratio T/R to be used in the Potvin decision schemes (power calculations 
in stage 1 and sample size estimation for stage 2).
}
  \item{CV}{
Coefficient of variation of the intra-subject variability as ratio.
}
  \item{targetpower}{
Power threshold in the power monitoring steps and 
power to achieve in the sample size estimation step.
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
  \item{usePE}{
If \code{TRUE} the sample size estimation step is done with PE and mse of 
stage 1.\cr
Defaults to \code{FALSE} in wich case the sample size is estimated with GMR and 
mse of stage 1 analogous to Potvin et. al.\cr
The power inspection steps in the Potvin methods are always done with the
GMR argument and mse (CV) of stage 1.
}
  \item{Nmax}{
Futility criterion. If set to a finite value all studies simulated in which a 
sample size >Nmax is obtained will be regarded as BE=FAIL.\cr
Set this argument to \code{Inf}, the default, to work without that futility
criterion.
}
  \item{min.n2}{
Minimum sample size of stage 2. Defaults to zero.\cr
If the sample size estimation step gives N <= n1 the sample size for stage 2 
will be set to \code{min.n2}, i.e. the total sample size to \code{n1+min.n2}.
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
Percentiles to be used for the presentation of the distribution of n(total)=n1+n2.\cr
Defaults to \code{c(0.05, 0.5, 0.95)} to obtain the 5\% and 95\% percentiles
and the median.
}
  \item{nsims}{
Number of studies to simulate.\cr 
If missing is set to 1E+05 = 100 000 or to 1E+06 = 1 mio if you are calculating 
'alpha', i.e. with \code{theta0} at border or outside acceptance range 
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
Potvin et al.\cr
The underlying subject data are assumed to be evaluated after log-transformation.
But instead of simulating subject data the statistics pe1, mse1 and pe2, SS2 are
simulated via their associated distributions (normal and chi-squared 
distri's).
}
\value{
Returns an object of class "pwrtsd" with all the input arguments and results 
as components.\cr
The class "pwrtsd" has an S3 print method.\cr
The results are in the components:
\item{pBE}{Contains the ratio of studies found BE.}
\item{pBE_s1}{Ratio of studies found BE in stage 1.}
\item{pct_s2}{Percentage of studies continuing to stage 2.}
\item{nmean}{Mean of n(total).}
\item{nrange}{Range (min, max) of n(total).}
\item{nperc}{Vector of percentiles of the distribution of n(total).}
\item{ntable}{Object of class "table" summarizing the discrete distribution of 
n(total) via its distinct values and counts of occurences of these values.\cr
This component is only given back if \code{usePE==FALSE} or otherwise if 
\code{is.finite(Nmax)}, i.e. a futility criterion is used.}
}
\references{
Potvin D et al.\cr
\emph{Sequential design approaches for bioequivalence studies with 
crossover designs}\cr
Pharm Stat 7(4), 245-62 (2008) \doi{10.1002/pst.294}

Montague TH et al.\cr
\emph{Additional results for 'Sequential design approaches for bio-equivalence studies 
with crossover designs'}\cr
Pharm Stat. 2012 Jan-Feb;11(1):8-13. Epub 2011 Feb 10. \doi{10.1002/pst.483}

Fuglsang A\cr
\emph{Controlling type I errors for two-stage bioequivalence study designs}\cr
Clinical Research and Regulatory Affairs 28(4), 100-105 (2011)\cr
\doi{10.3109/10601333.2011.631547}

Fuglsang A\cr
\emph{Sequential Bioequivalence Trial Designs with Increased Power and 
Controlled Type I Error Rates}\cr
AAPS J. 2013 Jul; 15(3):659-61. Epub 2013 Mar 30. \doi{10.1208/s12248-013-9475-5}

Fuglsang A\cr
\emph{Futility Rules in Bioequivalence Trials with Sequential Designs}\cr
AAPS J. 2014 Jan;16(1):79-82. Epub 2013 Nov 12. \doi{10.1208/s12248-013-9540-0}

\enc{Schütz}{Schuetz} H.\cr
\emph{Two-stage designs in bioequivalence trials}\cr
Eur J Clin Pharmacol 71(3):271-81 (2015)\cr
published online 22 Jan 2015 \doi{10.1007/s00228-015-1806-2}

Kieser M, Rauch G\cr
\emph{Two-stage designs for cross-over bioequivalence trials}\cr
Stat Med. 2015 Jul 20;34(16):2403-16\cr
Epub 2015 Mar 24 \doi{10.1002/sim.6487}

Zheng Ch, Zhao L, Wang J\cr
\emph{Modifications of sequential designs in bioequivalence trials}\cr
Pharm Stat 14 (3), 180-188, May/June 2015\cr
published online: 9 Feb 2015 \doi{10.1002/pst.1672}

Xu et al.\cr
\emph{Optimal adaptive sequential designs for crossover bioequivalence studies}\cr
Pharm Stat. 2016 Jan-Feb;15(1):15-27\cr
[Epub ahead of print] 2015 Nov 5 \doi{10.1002/pst.1721} 
}
\author{
D. Labes
}
\seealso{
\code{\link{power.2stage.p}} for analogous calculations for 2-group parallel 
design.\cr
\code{\link{power.2stage.fC}} for analogous calculations with futility check 
based on point estimate of stage 1.
}
\examples{
# using all the defaults and 24 subjects in stage 1, CV of 25\%
power.2stage(n1=24, CV=0.25)
# computation time ~ 1 sec
# 
# as above, but save results for further use
res <- power.2stage(n1=24, CV=0.25)
\dontrun{
# representation of the discrete distribution of n(total)
# via plot method of object with class "table" which creates a 'needle' plot
plot(res$ntable/sum(res$ntable), ylab="Density", main="Distribution of n(total)")
#
# If you prefer a histogram instead (IMHO not the preferred plot):
# reconstruct the ntotal values from the ntable
ntot <- rep.int(as.integer(names(res$ntable)), times=as.integer(res$ntable))
# annotated histogram
hist(ntot, freq=FALSE, breaks=res$nrange[2]-res$nrange[1])
abline(v=c(res$nmean, res$nperc[["50\%"]]), lty=c(1, 3))
legend("topright", box.lty=0, legend=c("mean", "median"), lty=c(1, 3), cex=0.9)}
}
