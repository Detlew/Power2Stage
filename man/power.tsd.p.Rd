\encoding{utf-8}
\name{power.tsd.p}
\alias{power.tsd.p}
\alias{power.tsd.pAF}
\alias{power.2stage.p}
\alias{power.2stage.pAF}
\title{
Power calculation of adaptive 2-stage BE studies in 2-group parallel designs
}
\description{
  This functions calculate the \sQuote{empirical} power of 2-stage BE studies
  with 2 parallel groups according to Fuglsang via simulations. The Potvin decision
  schemes are modified to include a futility criterion Nmax and to allow the
  sample size estimation step to be done with point estimate and variabilities
  from \ifelse{html}{\out{stage&nbsp;1}}{stage 1} (fully adaptive).\cr\cr
  Function \verb{power.tsd.pAF()} performes exactly as described in
  Fuglsang\enc{’}{'}s paper, namely the power monitoring steps and the sample
  size estimation are based always on the pooled \emph{t}-test.\cr
  Function \code{power.tsd.p()} with argument \verb{test="welch"} on the
  other hand uses the genuine power of Welch\enc{’}{'}s test. Moreover it
  accepts unequal \ifelse{html}{\out{stage&nbsp;1}}{stage 1} treatment groups.
  }
\usage{
power.tsd.p(method = c("B", "C"), alpha0 = 0.05, alpha = c(0.0294, 0.0294),
            n1, GMR, CV, targetpower = 0.8,
            pmethod = c("nct", "exact", "shifted"),
            usePE = FALSE, Nmax = Inf, test = c("welch", "t-test", "anova"),
            theta0, theta1, theta2, npct = c(0.05, 0.5, 0.95), nsims,
            setseed = TRUE, details = FALSE)

power.tsd.pAF(method = c("B", "C"), alpha0 = 0.05, alpha = c(0.0294, 0.0294),
              n1, GMR, CV, targetpower = 0.8,
              pmethod = c("shifted", "nct", "exact"),
              usePE = FALSE, Nmax = Inf, test = c("welch", "t-test", "anova"),
              theta0, theta1, theta2, npct = c(0.05, 0.5, 0.95), nsims,
              setseed = TRUE, details = FALSE)
}
\arguments{
  \item{method}{
    Decision schemes according to Potvin \emph{et.al.} (defaults to \verb{"B"}).\cr
    Potvin\enc{’}{'}s \sQuote{method D} can be obtained by choosing \verb{"C"} but
    setting \verb{alpha=c(0.028, 0.028)}.
  }
  \item{alpha0}{
    Alpha value for the first step(s) in Potvin \verb{"C"}, the power inspection
    and BE decision if power > targetpower.\cr
    Defaults to 0.05.
  }
  \item{alpha}{
    Vector (two elements) of the nominal alphas for the two stages.\cr
    Defaults to Pocock\enc{’}{'}s alpha setting \verb{alpha=c(0.0294, 0.0294)}.
  }
  \item{n1}{
    (Total) sample size of \ifelse{html}{\out{stage&nbsp;1}}{stage 1}.\cr
    Function \verb{power.tsd.p()} accepts also a vector of length 2 with
    the numbers in treatment groups of \ifelse{html}{\out{stage&nbsp;1}}{stage 1}.\cr
    If given as a scalar the total \verb{n1} should be even.
  }
  \item{GMR}{
    Ratio T/R to be used in decision scheme (power calculations in
    \ifelse{html}{\out{stage&nbsp;1}}{stage 1} and sample size estimation for
    \ifelse{html}{\out{stage&nbsp;2}}{stage 2}).
  }
  \item{CV}{
    Coefficient of variation of the total variability
    (use \emph{e.g.}, 0.3 for 30\%)\cr
    Can be a vector with two elements. In that case \verb{CV[1]} is for
    Test and \verb{CV[2]} for the Reference treatment.
  }
  \item{targetpower}{
    Power threshold in the power monitoring steps and power to achieve
    in the sample size estimation step.
  }
  \item{pmethod}{
    Power calculation method, also to be used in the sample size estimation for
    \ifelse{html}{\out{stage&nbsp;2}}{stage 2}.\cr
    Implemented are \verb{"nct"} (approximate calculations via non-central
    \emph{t}-distribution, \verb{"exact"} (exact calculations via Owen\enc{’}{'}s Q),
    and \verb{"shifted"} (approximate calculation via shifted central
    \emph{t}-distribution\cr
    Defaults to \verb{"nct"} for speed reasons in function \verb{power.tsd.p()}
    as a reasonable compromise between speed and accuracy in the sample size
    estimation step.\cr
    Defaults to \verb{"shifted"} in function \verb{power.tsd.pAF()} for
    consistancy with Fuglsang.\cr
  }
  \item{usePE}{
    If \verb{TRUE} the sample size estimation step is done with \acronym{MSE}
    \bold{and} \acronym{PE} of \ifelse{html}{\out{stage&nbsp;1}}{stage 1}.\cr
    Defaults to \code{FALSE} \emph{i.e.}, the sample size is estimated with
    \verb{GMR} and \acronym{MSE} (calculated from \verb{CV}) of
    \ifelse{html}{\out{stage&nbsp;1}}{stage 1} analogous to Potvin \emph{et. al}.\cr
    NB: The power inspection steps in the Potvin methods are always done with the
    \verb{GMR} argument and \acronym{MSE} (\verb{CV}) of
    \ifelse{html}{\out{stage&nbsp;1}}{stage 1}.
  }
  \item{Nmax}{
    Futility criterion. If set to a finite value, all studies simulated in which a
    sample size \verb{>Nmax} is obtained will be regarded as BE=FAIL.\cr
    Set this argument to \verb{Inf}, the default, to work without that futility
    criterion.
  }
  \item{test}{
    Test on which the CI calculations are based on.\cr
    Defaults to \verb{"welch"} = Welch\enc{’}{'}s \emph{t}-test accounting for
    heteroscedasticity in the variabilities of Test and Reference, but neglecting
    stage effects.\cr
    \verb{"anova"} calculates the
    \ifelse{html}{\out{100(1&ndash;2&alpha;)}}{100(1–2\eqn{\alpha}{alpha})}
    confidence interval based on an \acronym{ANOVA} with treatment and stage
    in the model.\cr
    \verb{"t-test"} calculates the
    \ifelse{html}{\out{100(1&ndash;2&alpha;)}}{100(1–2\eqn{\alpha}{alpha})}
    confidence interval based on the \emph{t}-test assuming equal
    variabilities of Test and Reference and neglecting stage effects.
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
  The calculations follow in principle the simulations as described by Fuglsang.\cr
  The underlying subject data are assumed to be evaluated after log-transformation.
  But instead of simulating subject data the statistics (mean and variance of Test
  and Reference of \ifelse{html}{\out{stage&nbsp;1}}{stage 1} and
  \ifelse{html}{\out{stage&nbsp;2}}{stage 2}) are simulated via their associated
  distributions (normal and \ifelse{html}{\out{&chi;<sup>2</sup>}}{\eqn{\chi^{2}}}).
}
\value{
  Returns an object of class \verb{"pwrtsd"} with all the input arguments and results
  as components.\cr
  The class \verb{"pwrtsd"} has an S3 print method.\cr
  The results are in the components:
  \item{pBE}{Fraction of studies found BE.}
  \item{pBE_s1}{Fraction of studies found BE in \ifelse{html}{\out{stage&nbsp;1}}{stage 1}.}
  \item{pct_s2}{Percentage of studies continuing to \ifelse{html}{\out{stage&nbsp;2}}{stage 2}.}
  \item{nmean}{Mean of n(total).}
  \item{nrange}{Range (min, max) of n(total).}
  \item{nperc}{Percentiles of the distribution of n(total).}
  \item{ntable}{Object of class \verb{"table"} summarizing the discrete distribution of
  n(total) via its distinct values and counts of occurences of these values.
  This component is only given back if \verb{usePE==FALSE} or otherwise if
  \verb{is.finite(Nmax)}, \emph{i.e.}, a futility criterion is used.}
}
\author{
D. Labes
}
\references{
Fuglsang A. \emph{Sequential Bioequivalence Approaches for Parallel Design.}\cr
AAPS J. 2014; 16(3):373--8. \doi{10.1208/s12248-014-9571-1}

Potvin D, DiLiberti CE, Hauck WW, Parr AF, Schuirmann DJ, Smith RA. \emph{Sequential design approaches for bioequivalence studies with crossover designs.}\cr
Pharm Stat. 2008; 7(4):245--62. \doi{10.1002/pst.294}
}
\seealso{
\code{\link{power.2stage}} for analogous calculations for the 2\enc{×}{x}2 crossover.
}
\examples{
# using all the defaults
power.tsd.p(n1=48, CV=0.25)
}
