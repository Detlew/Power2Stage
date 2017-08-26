\encoding{utf-8}
\name{power.2stage}
\alias{power.2stage}
\title{
Power calculation of adaptive 2-stage BE studies in a 2x2 crossover design
}
\description{
  This function calculates the \sQuote{empiric} power of 2-stage BE studies according
  to Potvin \emph{et al.} via simulations. The Potvin methods are modified to include
  a futility criterion Nmax and to allow the sample size estimation step to be
  done with the point estimate (\acronym{PE}) and \acronym{MSE} (calculated from CV)
  of \ifelse{html}{\out{stage&nbsp;1}}{stage~1}.
  }
\usage{
power.2stage(method = c("B", "C", "B0"), alpha0 = 0.05,
             alpha = c(0.0294, 0.0294), n1, GMR, CV, targetpower = 0.8,
             pmethod = c("nct", "exact", "shifted"),
             usePE = FALSE, Nmax = Inf, min.n2 = 0, theta0, theta1, theta2,
             npct = c(0.05, 0.5, 0.95), nsims, setseed = TRUE, details = FALSE)
}
\arguments{
  \item{method}{
    Decision schemes according to Potvin \emph{et.al.} (defaults to \verb{"B"}).\cr
    Potvin\enc{’}{'}s \sQuote{Method D} can be obtained by choosing \verb{"C"} but setting
    \verb{alpha=c(0.028, 0.028)}.\cr
    Montague\enc{’}{'}s \sQuote{Method D} can be obtained by choosing \verb{"C"} but setting
    \verb{alpha=c(0.028, 0.028)} and \verb{GMR=0.9}.\cr
    \verb{method="B0"} uses the decision scheme of Zheng \emph{et al.} \acronym{MSDBE}
    (modified sequential design for \acronym{B}E studies) which differs from B in case
    of different alphas w.r.t. power monitoring and BE decision in case of power
    >= target power.
  }
  \item{alpha0}{
    Alpha value for the first step(s) in Potvin \verb{"C"}, the power inspection
    and BE decision if power > targetpower. Defaults to 0.05.\cr
    Only observed if \verb{method="C"}.
  }
  \item{alpha}{
    Vector (two elements) of the nominal alphas for the two stages. Defaults to
    Pocock\enc{’}{'}s setting \verb{alpha=c(0.0294, 0.0294)}.\cr
    Common values together with other arguments are:\cr
    \verb{rep(0.0294, 2)}: Potvin \emph{et al.} \sQuote{Method B}\cr
    \verb{rep(0.0269, 2)}: Fulgsang \sQuote{Method C/D} \verb{(method="C", GMR=0.9, targetpower=0.9)}\cr
    \verb{rep(0.0274, 2)}: Fuglsang \sQuote{Method C/D} \verb{(method="C", targetpower=0.9)}\cr
    \verb{rep(0.0280, 2)}: Montague \emph{et al.} \sQuote{Method D} \verb{(method="C", GMR=0.9)}\cr
    \verb{rep(0.0284, 2)}: Fulgsang \sQuote{Method B} \verb{(GMR=0.9, targetpower=0.9)}\cr
    \verb{rep(0.0304, 2)}: Kieser & Rauch\cr
    \verb{c(0.01, 0.04)}: Zheng \emph{et al.} \sQuote{MSDBE} \verb{(method="B0")}
  }
  \item{n1}{
    Sample size of \ifelse{html}{\out{stage&nbsp;1}}{stage~1}.
  }
  \item{GMR}{
    Ratio T/R to be used in decision scheme (power calculations in
    \ifelse{html}{\out{stage&nbsp;1}}{stage~1} and sample size estimation for
    \ifelse{html}{\out{stage&nbsp;2}}{stage~2}).
  }
  \item{CV}{
    Coefficient of variation of the intra-subject variability
    (use \emph{e.g.}, 0.3 for 30\%).
  }
  \item{targetpower}{
    Power threshold in the power monitoring steps and power to achieve
    in the sample size estimation step.
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
  \item{usePE}{
    If \verb{TRUE} the sample size estimation step is done with \acronym{MSE}
    \bold{and} \acronym{PE} of \ifelse{html}{\out{stage&nbsp;1}}{stage~1}.\cr
    Defaults to \code{FALSE} \emph{i.e.}, the sample size is estimated with
    \verb{GMR} and \acronym{MSE} (calculated from \verb{CV}) of
    \ifelse{html}{\out{stage&nbsp;1}}{stage~1} analogous to Potvin \emph{et. al}.\cr
    NB: The power inspection steps in the Potvin methods are always done with the
    \verb{GMR} argument and \acronym{MSE} (\verb{CV}) of
    \ifelse{html}{\out{stage&nbsp;1}}{stage~1}.
  }
  \item{Nmax}{
    Futility criterion. If set to a finite value, all studies simulated in which a
    sample size \verb{>Nmax} is obtained will be regarded as BE=FAIL.\cr
    Set this argument to \verb{Inf}, the default, to work without that futility
    criterion.
  }
  \item{min.n2}{
    Minimum sample size of \ifelse{html}{\out{stage&nbsp;2}}{stage~2}. Defaults to zero.\cr
    If the sample size estimation step gives \verb{N < n1} the sample size for
    \ifelse{html}{\out{stage&nbsp;2}}{stage~2} will be forced to \verb{min.n2},
    \emph{i.e.}, the total sample size to \code{n1+min.n2}.
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
  The class \verb{"pwrtsd"}" has an S3 print method.\cr
  The results are in the components:
  \item{pBE}{Fraction of studies found BE.}
  \item{pBE_s1}{Fraction of studies found BE in \ifelse{html}{\out{stage&nbsp;1}}{stage~1}.}
  \item{pct_s2}{Percentage of studies continuing to \ifelse{html}{\out{stage&nbsp;2}}{stage~2}.}
  \item{nmean}{Mean of n(total), aka average total sample size (\acronym{ASN}).}
  \item{nrange}{Range (min, max) of n(total).}
  \item{nperc}{Vector of percentiles of the distribution of n(total).}
  \item{ntable}{Object of class \verb{"table"} summarizing the discrete distribution of
  n(total) via its distinct values and counts of occurences of these values.\cr
  This component is only given back if \verb{usePE==FALSE} or otherwise if
  \verb{is.finite(Nmax)}, \emph{i.e.}, a futility criterion is used.}
}
\author{
D. Labes
}
\references{
Potvin D, DiLiberti CE, Hauck WW, Parr AF, Schuirmann DJ, Smith RA. \emph{Sequential design approaches for bioequivalence studies with crossover designs.}\cr
Pharm Stat. 2008; 7(4):245--62. \doi{10.1002/pst.294}

Montague TH, Potvin D, DiLiberti CE, Hauck WW, Parr AF, Schuirmann DJ. \emph{Additional results for \sQuote{Sequential design approaches for bioequivalence studies with crossover designs}.}\cr
Pharm Stat. 2011; 11(1):8--13. \doi{10.1002/pst.483}

Fuglsang A. \emph{Controlling type I errors for two-stage bioequivalence study designs.}\cr
Clin Res Reg Aff. 2011; 28(4):100--5. \doi{10.3109/10601333.2011.631547}

Fuglsang A. \emph{Sequential Bioequivalence Trial Designs with Increased Power and
Controlled Type I Error Rates.}\cr
AAPS J. 2013; 15(3):659--61. \doi{10.1208/s12248-013-9475-5}

Fuglsang A. \emph{Futility Rules in Bioequivalence Trials with Sequential Designs.}\cr
AAPS J. 2014; 16(1):79--82. \doi{10.1208/s12248-013-9540-0}

\enc{Schütz}{Schuetz} H. \emph{Two-stage designs in bioequivalence trials.}\cr
Eur J Clin Pharmacol. 2015; 71(3):271--81. \doi{10.1007/s00228-015-1806-2}

Kieser M, Rauch G. \emph{Two-stage designs for cross-over bioequivalence trials.}\cr
Stat Med. 2015; 34(16):2403--16. \doi{10.1002/sim.6487}

Zheng Ch, Zhao L, Wang J. \emph{Modifications of sequential designs in bioequivalence trials.}\cr
Pharm Stat. 2015; 14(3):180--8. \doi{10.1002/pst.1672}

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
# via plot method of object with class "table" which creates a
# 'needle' plot
plot(res$ntable/sum(res$ntable), ylab="Density",
     xlab=expression("n"[total]), las=1,
     main=expression("Distribution of n"[total]))
#
# If you prefer a histogram instead (IMHO, not the preferred plot):
# reconstruct the ntotal values from the ntable
ntot <- rep.int(as.integer(names(res$ntable)),
                times=as.integer(res$ntable))
# annotated histogram
hist(ntot, freq=FALSE, breaks=res$nrange[2]-res$nrange[1],
     xlab=expression("n"[total]), las=1,
     main=expression("Histogram of n"[total]))
abline(v=c(res$nmean, res$nperc[["50\%"]]), lty=c(1, 3))
legend("topright", box.lty=0, legend=c("mean", "median"),
       lty=c(1, 3), cex=0.9)}
}
