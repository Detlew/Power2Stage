\encoding{utf-8}
\name{power.tsd.GS}
\alias{power.tsd.GS}
\alias{power.2stage.GS}
\title{
Power calculation for non-adaptive group sequential (2-stage) BE studies
}
\description{
  This function calculates the \sQuote{empiric} power of group sequential
  2-stage \acronym{BE} in 2\enc{×}{x}2 crossover designs via simulations.\cr
  The number of subjects in both stages has to be prespecified (non-adaptive).
  }
\usage{
power.tsd.GS(alpha = c(0.0294, 0.0294), n, CV, theta0, theta1, theta2,
                fCrit = c("CI", "PE"), fClower, fCupper, nsims,
                setseed = TRUE, details = FALSE)
}
\arguments{
  \item{alpha}{
    Vector of the two nominal alpha values to be used in the
    \ifelse{html}{\out{100(1&ndash;2&alpha;)}}{100(1–2\eqn{\alpha}{alpha})}
    confidence interval calculations in the two stages.\cr
    Use something like package \verb{ldbounds} for choosing the nominal alphas.
  }
  \item{n}{
    Vector of the two sample sizes in \ifelse{html}{\out{stage&nbsp;1}}{stage~1}
    and \ifelse{html}{\out{stage&nbsp;2}}{stage~2}. \verb{n(total)} is \verb{n[1]+n[2]}
    if a second stage is necessary. Otherwise it is \verb{n[1]}.
  }
  \item{CV}{
    Coefficient of variation of the intra-subject variability
    (use \emph{e.g.}, 0.3 for 30\%).
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
  \item{fCrit}{
    Futility criterion.\cr
    If set to \verb{"PE"} the study stops after
    \ifelse{html}{\out{stage&nbsp;1}}{stage~1} if not \acronym{BE} and if the point
    estimate (\acronym{PE}) of \ifelse{html}{\out{stage&nbsp;1}}{stage~1} evaluation
    is outside the range defined in the next two  arguments \verb{"fClower"} and
    \verb{"fCupper"}.\cr
    If set to \verb{"CI"} the study stops after \ifelse{html}{\out{stage&nbsp;1}}{stage~1}
    if not \acronym{BE} and if the 90\% confidence interval of
    \ifelse{html}{\out{stage&nbsp;1}}{stage~1} evaluation is outside the range
    defined in the next two arguments.\cr
    Defaults to \verb{"CI"}.
  }
  \item{fClower}{
    Lower limit of the futility criterion. Defaults to \verb{0} if missing,
    \emph{i.e.}, no futility criterion.
  }
  \item{fCupper}{
    Upper limit of the futility criterion. Defaults to \verb{1/fClower} if missing.
  }
  \item{nsims}{
    Number of studies to simulate.\cr
    If missing, \verb{nsims} is set to 1E+05 = 100,000 or to 1E+06 = 1 Mio if
    estimating the empiric Type I Error (\verb{'alpha'}), \emph{i.e.}, with \verb{theta0}
    at the border or outside the acceptance range \verb{theta1} \ldots \verb{theta2}.
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
  \emph{et al.} for adaptive designs, but with no interim power inspection and
  pre-specified (fixed) number of subjects in \ifelse{html}{\out{stage&nbsp;2}}{stage~2}.\cr
  The underlying subject data are assumed to be evaluated after log-transformation.
  But instead of simulating subject data, the statistics pe1, mse1 and pe2, SS2 are
  simulated via their associated distributions (normal and
  \ifelse{html}{\out{&chi;<sup>2</sup>}}{\eqn{\chi^{2}}} distributions).
  }
\value{
  Returns an object of class \verb{"pwrtsd"} with all the input arguments and results
  as components.\cr
  The class \verb{"pwrtsd"} has a S3 print method.\cr
  The results are in the components:
  \item{pBE}{Fraction of studies found BE.}
  \item{pBE_s1}{Fraction of studies found BE in \ifelse{html}{\out{stage&nbsp;1}}{stage~1}.}
  \item{pct_s2}{Percentage of studies continuing to \ifelse{html}{\out{stage&nbsp;2}}{stage~2}.}
  }
\note{
The code is reasonable fast. 1E6 sims take ~ 1 sec on my machine.
Even 1E7 sims are meanwhile possible without too much beer. :-)\cr\cr
}
\author{
D. Labes
}
\references{
Gould AL. \emph{Group sequential extensions of a standard bioequivalence testing procedure.}\cr
J Pharmacokin Biopharm. 1995; 23(1):57--86 \doi{10.1007/BF02353786}

Patterson SD, Jones B. \emph{Bioequivalence and Statistics in Clinical Pharmacology.}\cr
Boca Raton: CRC Press; 2\ifelse{html}{\out{<sup>nd</sup>}}{\eqn{^{nd}}} edition 2016. Chapter 5.6 Optional Designs.

\enc{Schütz}{Schuetz} H. \emph{Two-stage designs in bioequivalence trials.}\cr
Eur J Clin Pharmacol. 2015; 71(3):271--81. \doi{10.1007/s00228-015-1806-2}

Kieser M, Rauch G. \emph{Two-stage designs for cross-over bioequivalence trials.}\cr
Stat Med. 2015; 34(16):2403--16. \doi{10.1002/sim.6487}

Zheng Ch, Zhao L, Wang J. \emph{Modifications of sequential designs in bioequivalence trials.}\cr
Pharm Stat. 2015; 14(3):180--8. \doi{10.1002/pst.1672}
}
\seealso{
\code{\link{power.tsd}} and \code{\link{power.tsd.p}} for adaptive sequential designs.
}
\examples{
# using the Pocock alpha settings and no futility rule
# (defaults), a CV of 20\% and 12 subjects in both stages,
# midway interim
power.tsd.GS(CV=0.2, n=c(12,12))
}
