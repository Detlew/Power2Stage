\encoding{utf-8}
\name{power.2stage.GS}
\alias{power.2stage.GS}
\title{
Power calculation for non-adaptive group sequential (2-stage) BE studies
}
\description{
This function calculates the 'empirical' power of group sequential 2-stage BE 
with 2x2 crossover designs via simulations.\cr
The number of subjects in both stages have to be prespecified (non-adaptive).
}
\usage{
power.2stage.GS(alpha = c(0.0294, 0.0294), n, CV, theta0, theta1, theta2, 
                fCrit = c("PE", "CI"), fClower, fCupper, nsims, setseed = TRUE, 
                details = FALSE)
}
\arguments{
  \item{alpha}{
Vector of the two nominal alpha values to be used in the 1-2*alpha confidence
intervals calculations for the two stages.\cr
Use something like package \code{ldbounds} for choosing the nominal alpha's.
}
  \item{n}{
Vector of the two sample sizes in stage 1 and stage 2.\cr
n(total) is n[1]+n[2] if a second stage is necessary. Else it is n[1].
}
  \item{CV}{
Coefficient of variation of the intra-subject variability as ratio.
}
  \item{theta0}{
True ratio of T/R for simulating. Defaults to 0.95 if missing.
}
  \item{theta1}{
Lower bioequivalence limit. Defaults to 0.8.
}
  \item{theta2}{
Upper bioequivalence limit. Defaults to 1.25.
}
  \item{fCrit}{
Futility criterion.\cr
If set to \code{"PE"} the study stops after stage 1 if not BE and if the point 
estimator of stage 1 evaluation is outside the range defined in the next two 
argumenmts \code{"fClower"} and \code{"fCupper"}.\cr
If set to \code{"CI"} the study stops after stage 1 if not BE and if the 
confidence interval of stage 1 evaluation is outside the range defined in the 
next two argumenmts.\cr
Defaults to \code{"PE"}.
}
  \item{fClower}{
Lower limit of the futility criterion.\cr
Defaults to \code{0} if missing, i.e. no futility criterion.
}
  \item{fCupper}{
Upper limit of the futility criterion.\cr
Defaults to \code{1/fClower} if missing.
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
The calculations follow in principle the simulations as described in 
Potvin et al. for adaptive designs, but with no interim power inspection and 
pre-specified (fixed) number of subjects in stage 2.\cr
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
\item{pBE_s1}{Ratio of studies found BE in stage 1.}
\item{pct_s2}{Percentage of studies continuing to stage 2.}
}
\references{
Gould, A.L. \cr
\emph{Group sequential extensions of a standard bioequivalence testing procedure}\cr
J Pharmacokinet Biopharm. 1995 Feb;23(1):57-86 \doi{10.1007/BF02353786}

Patterson S and Jones B\cr
\emph{BIOEQUIVALENCE and STATISTICS in CLINICAL PHARMACOLOGY}\cr
Chapman & Hall/CRC\cr
Boca Raton 2006\cr
Chapter 5.6 Optional Designs

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

}
\author{
D. Labes
}
\note{
The code is reasonable fast. 1E6 sims take ~ 1 sec on my machine.
Even 1E7 sims are possible without meanwhile too much beer :-)).\cr\cr
}
\seealso{
\code{\link{power.2stage}} and \code{\link{power.2stage.p}} for adaptive sequential designs.
}
\examples{
# using the Pocock alpha settings and no futility rule
# (defaults), a CV of 20\% and 12 subjects in both stages,
# midway interim 
power.2stage.GS(CV=0.2, n=c(12,12))
}
