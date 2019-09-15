\encoding{utf-8}
\name{power.tsd.fC}
\alias{power.tsd.fC}
\alias{power.2stage.fC}
\title{
Power calculation of adaptive 2-stage BE studies (2x2 crossover) with a
futility criterion for the point estimate of T/R or its 90\% CI
}
\description{
  This function calculates the \sQuote{empiric} power of 2-stage BE studies according
  to Potvin \emph{et al.} \sQuote{method B/C} via simulations. The Potvin methods are
  modified to include a futility criterion for the point estimate or for its
  90\%CI and to allow the sample size estimation step to be done with the point
  estimate (\acronym{PE}) and \acronym{MSE} of \ifelse{html}{\out{stage&nbsp;1}}{stage 1}.
  }
\usage{
power.tsd.fC(method = c("B", "C", "B0"), alpha0 = 0.05, alpha = c(0.0294, 0.0294),
             n1, CV, GMR, targetpower = 0.8, pmethod = c("nct", "exact", "shifted"),
             usePE = FALSE, powerstep = TRUE, min.n2=0, max.n=Inf,
             fCrit=c("CI", "PE"), fClower, fCupper, theta0, theta1, theta2,
             npct = c(0.05, 0.5, 0.95), nsims, setseed = TRUE, details = FALSE)
}
\arguments{
  \item{method}{
    Decision schemes according to Potvin \emph{et.al.} (defaults to \verb{"B"}).\cr
    Montague\enc{’}{'}s \sQuote{Method D} can be obtained by choosing \verb{"C"} but setting
    \verb{alpha=c(0.028, 0.028)}.\cr
    \sQuote{Method E} of Xu \emph{et al.} can be obtained by choosing \verb{"B"} and setting alphas,
    futility criterion \verb{"CI"}, \verb{max.n}, and \verb{n1} according to the reference.\cr
    \sQuote{Method F} can be obtained choosing \verb{"C"} with the appropriate design setting
    according to the reference.\cr
    \verb{method="B0"} uses the decision scheme of Zheng \emph{et al.} \acronym{MSDBE}
    (modified sequential design for \acronym{B}E studies) which differs from B in case
    of different alphas w.r.t. power monitoring and BE decision in case of power
    >= target power.
  }
  \item{alpha0}{
    Alpha value for the first step(s) in Potvin \verb{"C"}, the power inspection
    and BE decision if power > targetpower. Defaults to 0.05.\cr
    Only observed if \verb{method="C"}
  }
  \item{alpha}{
    Vector (two elements) of the nominal alphas for the two stages. Defaults to
    Pocock\enc{’}{'}s setting \verb{alpha=c(0.0294, 0.0294)}.\cr
    Common values together with other arguments are:\cr
    \verb{rep(0.0294, 2)}: Potvin \emph{et al.} \sQuote{Method B} \verb{(fCrit="CI", fCupper=Inf)}\cr
    \verb{rep(0.0269, 2)}: Fulgsang \sQuote{Method C/D} \verb{(method="C", GMR=0.9, targetpower=0.9, fCrit="CI", fCupper=Inf)}\cr
    \verb{rep(0.0274, 2)}: Fuglsang \sQuote{Method C/D} \verb{(method="C", targetpower=0.9, fCrit="CI", fCupper=Inf)}\cr
    \verb{rep(0.0280, 2)}: Montague \emph{et al.} \sQuote{Method D} \verb{(method="C", GMR=0.9, fCrit="CI", fCupper=Inf)}\cr
    \verb{rep(0.0284, 2)}: Fulgsang \sQuote{Method B} \verb{(GMR=0.9, targetpower=0.9, fCrit="CI", fCupper=Inf)}\cr
    \verb{rep(0.0304, 2)}: Kieser & Rauch \verb{(fCrit="CI", fCupper=Inf)}\cr
    \verb{c(0.01, 0.04)}: Zheng \emph{et al.} \sQuote{MSDBE} \verb{(method="B0", fCrit="CI", fCupper=Inf)}\cr
    \verb{c(0.0249, 0.0357)}: Xu \emph{et al.} \sQuote{Method E} for CV 10--30\% \verb{(fCrit="CI", fClower=0.9374, max.n=42)}\cr
    \verb{c(0.0254, 0.0363)}: Xu \emph{et al.} \sQuote{Method E} for CV 30--55\% \verb{(fCrit="CI", fClower=0.9305, max.n=42)}\cr
    \verb{c(0.0248, 0.0364)}: Xu \emph{et al.} \sQuote{Method F} for CV 10--30\% \verb{(method="C", fCrit="CI", fClower=0.9492, max.n=180)}\cr
    \verb{c(0.0259, 0.0349)}: Xu \emph{et al.} \sQuote{Method F} for CV 30--55\% \verb{(method="C", fCrit="CI", fClower=0.9350, max.n=180)}
  }
  \item{n1}{
    Sample size of \ifelse{html}{\out{stage&nbsp;1}}{stage 1}. For Xu\enc{’}{'}s
    methods the recommended sample size should be at least 18 (if CV 10--30\%) or
    48 (if CV 30--55\%).
  }
  \item{CV}{
    Coefficient of variation of the intra-subject variability
    (use \emph{e.g.}, 0.3 for 30\%).
  }
  \item{GMR}{
    Ratio T/R to be used in decision scheme (power calculations in
    \ifelse{html}{\out{stage&nbsp;1}}{stage 1} and sample size estimation for
    \ifelse{html}{\out{stage&nbsp;2}}{stage 2}).
  }
  \item{targetpower}{
    Power threshold in the power monitoring steps and
    power to achieve in the sample size estimation step.
  }
  \item{pmethod}{
    Power calculation method, also to be used in the sample size estimation for
    \ifelse{html}{\out{stage&nbsp;2}}{stage 2}.\cr
    Implemented are \verb{"nct"} (approximate calculations via non-central
    \emph{t}-distribution, \verb{"exact"} (exact calculations via Owen\enc{’}{'}s Q),
    and \verb{"shifted"} (approximate calculation via shifted central \emph{t}-distribution
    like in the paper of Potvin \emph{et al.}\cr
    Defaults to \verb{"nct"} as a reasonable compromise between speed and
    accuracy in the sample size estimation step.
  }
  \item{usePE}{
    If \verb{TRUE} the sample size estimation step is done with \acronym{MSE}
    \bold{and} PE of \ifelse{html}{\out{stage&nbsp;1}}{stage 1}.\cr
    Defaults to \verb{FALSE}, \emph{i.e.}, the sample size is estimated with anticipated
    (fixed) \verb{GMR} given as argument and \acronym{MSE} of
    \ifelse{html}{\out{stage&nbsp;1}}{stage 1} (analogous to Potvin \emph{et. al.}).
  }
  \item{powerstep}{
    If \verb{TRUE} (the default) the interim power monitoring step in
    \ifelse{html}{\out{stage&nbsp;1}}{stage 1} evaluation of \sQuote{method B}
    will be done as described in Potvin \emph{et.al.}\cr
    Setting this argument to \verb{FALSE} will omit this step.\cr
    Has no effect if \verb{method="C"} is choosen.
  }
  \item{min.n2}{
    Minimum sample size of \ifelse{html}{\out{stage&nbsp;2}}{stage 2}. Defaults to zero.\cr
    If the sample size estimation step gives \verb{N < n1} the sample size for
    \ifelse{html}{\out{stage&nbsp;2}}{stage 2} will be forced to \verb{min.n2},
    \emph{i.e.}, the total sample size to \code{n1+min.n2}.
  }
  \item{max.n}{
    If \verb{max.n} is set to a finite value the re-estimated total sample size (N) is
    set to \code{min(max.n, N)}.\cr
    Defaults to \verb{Inf} which is equivalent to not constrain the re-estimated sample size.\cr
    Attention! \verb{max.n} here is \bold{not} a futility criterion like \verb{Nmax}
    in other functions of this package.
  }
  \item{fCrit}{
    Futility criterion.\cr
    If set to \verb{"PE"} the study stops after
    \ifelse{html}{\out{stage&nbsp;1}}{stage 1} if not \acronym{BE} and if the point
    estimate (\acronym{PE}) of \ifelse{html}{\out{stage&nbsp;1}}{stage 1} evaluation
    is outside the range defined in the next two  arguments \verb{"fClower"} and
    \verb{"fCupper"}.\cr
    If set to \verb{"CI"} the study stops after \ifelse{html}{\out{stage&nbsp;1}}{stage 1}
    if not \acronym{BE} and if the confidence interval of
    \ifelse{html}{\out{stage&nbsp;1}}{stage 1} evaluation is outside the range
    defined in the next two arguments.\cr
    Defaults to \verb{"PE"}.
    Futility criterion to use for \verb{PE} or \verb{CI}.
  }
  \item{fClower}{
    Lower futility limit for the \verb{PE} or \verb{CI} of
    \ifelse{html}{\out{stage&nbsp;1}}{stage 1}.\cr
    If the \verb{PE} or \verb{CI} is outside \verb{fClower} \ldots \verb{fCupper}
    the study is stopped in the interim with the result FAIL (not BE).\cr
    May be missing. Defaults then to 0.8 if \verb{fCrit="PE"} or 0.925 if
    \verb{fCrit="CI"}.
  }
  \item{fCupper}{
    Upper futility limit for the \verb{PE} or \verb{CI} of
    \ifelse{html}{\out{stage&nbsp;1}}{stage 1}.\cr
    Will be set to \verb{1/fClower} if missing.
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
    of the simulation steps. Defaults to \verb{FALSE}.
  }
}
\details{
  The calculations follow in principle the simulations as described in Potvin
  \emph{et al}.\cr
  The underlying subject data are assumed to be evaluated after log-transformation.
  But instead of simulating subject data, the statistics pe1, mse1 and pe2, SS2 are
  simulated via their associated distributions (normal and
  \ifelse{html}{\out{&chi;<sup>2</sup>}}{\eqn{\chi^{2}}} distributions).
}
\value{
  Returns an object of class \verb{"pwrtsd"} with all the input arguments and results
  as components.\cr
  The class \verb{"pwrtsd"} has an S3 print method.\cr
  The results are in the components:
  \item{pBE}{Fraction of studies found BE.}
  \item{pBE_s1}{Fraction of studies found BE in \ifelse{html}{\out{stage&nbsp;1}}{stage 1}.}
  \item{pct_s2}{Percentage of studies continuing to \ifelse{html}{\out{stage&nbsp;2}}{stage 2}.}
  \item{nmean}{Mean of n(total), aka average total sample size (\acronym{ASN}).}
  \item{nrange}{Range (min, max) of n(total).}
  \item{nperc}{Percentiles of the distribution of n(total).}
  \item{ntable}{Object of class \verb{"table"} summarizing the discrete distribution of
  n(total) via its distinct values and counts of occurences of these values.\cr
  This component is only given back if \code{usePE==FALSE} or
  \code{usePE==TRUE & fClower>0 & is.finite(fCupper)}, \emph{i.e.}, a futility range is used.}
}
\author{
D. Labes
}
\references{
Potvin D, DiLiberti CE, Hauck WW, Parr AF, Schuirmann DJ, Smith RA. \emph{Sequential design approaches for bioequivalence studies with crossover designs.}\cr
Pharm Stat. 2008; 7(4):245--62. \doi{10.1002/pst.294}

Montague TH, Potvin D, DiLiberti CE, Hauck WW, Parr AF, Schuirmann DJ. \emph{Additional results for \sQuote{Sequential design approaches for bioequivalence studies with crossover designs}.}\cr
Pharm Stat. 2011; 11(1):8--13. \doi{10.1002/pst.483}

Fuglsang A. \emph{Sequential Bioequivalence Trial Designs with Increased Power and
Controlled Type I Error Rates.}\cr
AAPS J. 2013; 15(3):659--61. \doi{10.1208/s12248-013-9475-5}

\enc{Schütz}{Schuetz} H. \emph{Two-stage designs in bioequivalence trials.}\cr
Eur J Clin Pharmacol. 2015; 71(3):271--81. \doi{10.1007/s00228-015-1806-2}

Kieser M, Rauch G. \emph{Two-stage designs for cross-over bioequivalence trials.}\cr
Stat Med. 2015; 34(16):2403--16. \doi{10.1002/sim.6487}

Zheng Ch, Zhao L, Wang J. \emph{Modifications of sequential designs in bioequivalence trials.}\cr
Pharm Stat. 2015; 14(3):180--8. \doi{10.1002/pst.1672}

Xu J, Audet C, DiLiberti CE, Hauck WW, Montague TH, Parr TH, Potvin D, Schuirmann DJ.
\emph{Optimal adaptive sequential designs for crossover bioequivalence studies.}\cr
Pharm Stat. 2016;15(1):15--27. \doi{10.1002/pst.1721}
}
\seealso{
\code{\link{power.tsd}}
}
\examples{
# using all the defaults
power.tsd.fC(CV=0.25, n1=24)
# run-time ~1 sec
\dontrun{
# as above but storing the results
res <- power.tsd.fC(CV=0.25, n1=24)
# representation of the discrete distribution of n(total)
# via plot method of object with class "table" which creates a
# 'needle' plot
plot(res$ntable/sum(res$ntable), ylab="Density",
     xlab=expression("n"[total]), las=1,
     main=expression("Distribution of n"[total]))}
}
