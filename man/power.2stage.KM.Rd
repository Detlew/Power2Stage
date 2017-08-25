\encoding{utf-8}
\name{power.2stage.KM}
\alias{power.2stage.KM}
\title{
Power calculation of adaptive 2-stage BE studies (2x2 crossover)
(variant of power.2stage to obtain the results of Karalis / Macheras)
}
\description{
  This function calculates the empirical power of 2-stage BE studies according
  to Potvin \emph{et al.} via simulations. The Potvin methods are modified as
  described by Karalis & Macheras to include a futility criterion Nmax and to
  perform the power calculation steps and the sample size estimation step in
  the decision schemes with the \acronym{MSE} (calculated from CV) \strong{and}
  the point estimate (\acronym{PE}) of T/R from \ifelse{html}{\out{stage&nbsp;1}}{stage~1}.
  }
\usage{
power.2stage.KM(method = c("C", "B"), alpha0 = 0.05, alpha = c(0.0294, 0.0294),
                n1, CV, targetpower = 0.8, pmethod = c("nct", "exact"),
                Nmax = 150, theta0, theta1, theta2, npct = c(0.05, 0.5, 0.95),
                nsims, setseed = TRUE, details = FALSE)
}
\arguments{
  \item{method}{
    Decision schemes according to Potvin \emph{et al.}\cr
    Default is \verb{"C"} aka \acronym{TSD} in the paper of Karalis & Macheras
    if setting \verb{alpha=c(0.0294, 0.0294)}.\cr
    \acronym{TSD-1} of Karalis can be obtained by choosing \verb{"C"} but
    setting \verb{alpha=c(0.028, 0.028)}.\cr
    \acronym{TSD-2} of Karalis can be obtained by choosing \verb{"B"} and
    setting \verb{alpha=c(0.0294, 0.0294)}.\cr
  }
  \item{alpha0}{
    Alpha value for the first step(s) in Potvin C aka \acronym{TSD} of Karalis &
    Macheras or \acronym{TSD-1} of Karalis, the power inspection and BE decision
    if power > targetpower.\cr
    Defaults to 0.05.
  }
  \item{alpha}{
    Vector (two elements) of the nominal alphas for the two stages.\cr
    Defaults to Pocock\enc{’}{'}s alpha setting \verb{alpha=c(0.0294, 0.0294)}
    as in \acronym{TSD} of Karalis & Macheras.
  }
  \item{n1}{
    Sample size of \ifelse{html}{\out{stage&nbsp;1}}{stage~1}.
  }
  \item{CV}{
    Coefficient of variation of the intra-subject variability
    (use \emph{e.g.}, 0.3 for 30\%).
  }
  \item{targetpower}{
    Power threshold in the first step of Potvin \verb{"C"} and power to
    achieve in the sample size estimation step.
  }
  \item{pmethod}{
    Power calculation method, also to be used in the sample size estimation for
    \ifelse{html}{\out{stage&nbsp;2}}{stage~2}.\cr
    Implemented are "\verb{"nct"} (approximate calculations via non-central
    \emph{t}-distribution and \verb{"exact"} (exact calculations via Owen\enc{’}{'}s Q).\cr
    Defaults to \verb{"nct"} as a reasonable compromise between speed and
    accuracy in the sample size estimation step.
  }
  \item{Nmax}{
    Futility criterion. If set to a finite value all studies simulated in which
    a sample size >Nmax is obtained will be regarded as BE=FAIL. Defaults to 150,
    as recommended by Karalis & Macheras.\cr
    Set this argument to \verb{Inf}, to work without that futility criterion.
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
    Percentiles to be used for the presentation of the distribution of
    \verb{n(total)=n1+n2}.\cr
    Defaults to \verb{c(0.05, 0.5, 0.95)} to obtain the 5\% and 95\% percentiles
    and the median.
  }
  \item{nsims}{
    Number of studies to simulate.\cr
    If missing, \verb{nsims} is set to 1E+05 = 100,000 or to 1E+06 = 1 Mio
    if estimating the empiric Type I Error (\verb{'alpha'}), \emph{i.e.},
    with \verb{theta0} at the border or outside the acceptance range
    \verb{theta1} \ldots \verb{theta2}.
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
  \ifelse{html}{\out{&chi;<sup>2</sup>}}{\eqn{\chi^{2}}} distributions).\cr\cr
  In contrast to Potvin \emph{et al.} the power calculation steps as well as the
  sample size adaption step of the decision schemes are done using the \acronym{MSE}
  (calculated from CV) \bold{and} the point estimate from
  \ifelse{html}{\out{stage&nbsp;1}}{stage~1}.\cr
  This resembles the methods described in Karalis & Macheras and Karalis.
  }
\value{
  Returns an object of class \verb{"pwrtsd"} with all the input arguments and results
  as components.\cr
  The class \verb{"pwrtsd"} has a S3 print method.\cr
  The results are in the components:
  \item{pBE}{Fraction of studies found BE.}
  \item{pBE_s1}{Fraction of studies found BE in \ifelse{html}{\out{stage&nbsp;1}}{stage~1}.}
  \item{pct_s2}{Percentage of studies continuing to \ifelse{html}{\out{stage&nbsp;2}}{stage~2}.}
  \item{nmean}{Mean of n(total).}
  \item{nrange}{Range (min, max) of n(total).}
  \item{nperc}{Percentiles of the distribution of n(total).}
  \item{ntable}{Object of class \verb{"table"} summarizing the discrete distribution
  of n(total) via its distinct values and counts of occurences of these values.\cr
  This component is only given back if \verb{is.finite(Nmax)}.}
}
\author{
D. Labes
}
\references{
Potvin D, DiLiberti CE, Hauck WW, Parr AF, Schuirmann DJ, Smith RA. \emph{Sequential design approaches for bioequivalence studies with crossover designs.}\cr
Pharm Stat. 2008; 7(4):245--62. \doi{10.1002/pst.294}

Karalis V, Macheras P. \emph{An Insight into the Properties of a Two-Stage Design in Bioequivalence Studies.}\cr
Pharm Res. 2013; 30(7):1824--35. \doi{10.1007/s11095-013-1026-3}

Karalis V. \emph{The role of the upper sample size limit in two-stage bioequivalence designs.}\cr
Int J Pharm. 2013; 456(1):87--94. \doi{10.1016/j.ijpharm.2013.08.013}

Fuglsang A. \emph{Futility Rules in Bioequivalence Trials with Sequential Designs.}\cr
AAPS J. 2014; 16(1):79--82. \doi{10.1208/s12248-013-9540-0}

\enc{Schütz}{Schuetz} H. \emph{Two-stage designs in bioequivalence trials.}\cr
Eur J Clin Pharmacol. 2015; 71(3):271--81. \doi{10.1007/s00228-015-1806-2}
}
\examples{
# using all the defaults
# but too low number of sims to complain with the CRAN policy:
# "check time only a few seconds per example"
# minimum number of sims should be 1E5 for power, 1E6 sims for 'alpha'
power.2stage.KM(n1=16, CV=0.2, nsims=1E4)
# ~3 sec if nsims=1E5
}
