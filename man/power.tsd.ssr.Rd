\encoding{utf-8}
\name{power.tsd.ssr}
\alias{power.2stage.ssr}
\alias{power.tsd.ssr}
\title{
Power of 2-stage BE studies in 2x2 crossover designs with interim
sample size re-estimation
}
\description{
  This function calculates the \sQuote{empiric} power (via simulations) of
  2-stage BE studies with interim sample size re-estimation (\emph{i.e.},
  but no BE decision after \ifelse{html}{\out{stage&nbsp;1}}{stage~1}).
  The sample size re-estimation can be done blinded or unblinded.
  }
\usage{
power.tsd.ssr(alpha = 0.05, n1, GMR, CV, targetpower = 0.8,
              pmethod = c("nct", "exact", "shifted", "ls"), blind = FALSE,
              usePE = FALSE, min.n = 0, max.n = Inf, theta0, theta1, theta2,
              npct = c(0.05, 0.5, 0.95), nsims, setseed = TRUE, details = FALSE)
}
\arguments{
  \item{alpha}{
    Nominal type I error. Has to be adjusted in case of inflation of
    the Type I Error.
  }
  \item{n1}{
    Sample size of \ifelse{html}{\out{stage&nbsp;1}}{stage~1}.
  }
  \item{GMR}{
    Ratio T/R to be used in the sample size re-estimation.\cr
    Defaults to 0.95 if missing.
  }
  \item{CV}{
    Coefficient of variation of the intra-subject variability
    (use \emph{e.g.}, 0.3 for 30\%).\cr
    Anticipated population value.
  }
  \item{targetpower}{
    Power to achieve in the sample size estimation step.
  }
  \item{pmethod}{
    Power calculation method to be used in the sample size re-estimation for
    \ifelse{html}{\out{stage&nbsp;2}}{stage~2}.\cr
    Implemented are \verb{"nct"} (approximate calculations via non-central
    \emph{t}-distribution, \verb{"exact"} (exact calculations via Owen\enc{’}{'}s
    Q), and \verb{"shifted"} (approximate calculation via shifted central
    \emph{t}-distribution.\cr
    Also implemented is the large sample approximation as used in the references.\cr
    Defaults to \verb{"nct"} as a reasonable compromise between speed and
    accuracy in the sample size estimation step.
  }
  \item{blind}{
    If \verb{TRUE} the blinded estimate of the intra-subject variance, \emph{i.e.},
    the estimate from the period differences, is used in sample size estimation.\cr
    If \verb{FALSE} the usual \acronym{MSE} from \ifelse{html}{\out{stage&nbsp;1}}{stage~1}
    is used.\cr
    Defaults to \verb{FALSE} since most \acronym{BE} studies are open.
  }
  \item{usePE}{
    If \verb{TRUE} the point estimate from the interim analysis is used in the
    sample size re-estimation. Defaults to \verb{FALSE}.\cr
    \verb{usePE = TRUE} doesn\enc{’}{'}t make sense if \verb{blind = TRUE}. In that
    case the function issues a warning and \verb{usePE} is reset to \verb{usePE = FALSE}.
  }
  \item{min.n}{
    If \verb{min.n>n1}, the re-estimated sample size (\verb{N}) is set to
    \verb{max(min.n,N)}.\cr
    If \verb{min.n=0} (the default), no minimal sample size is applied.
  }
  \item{max.n}{
    If \verb{max.n} is set to a finite value the re-estimated sample size (\verb{N})
    is set to \verb{min(max.n,N)}.\cr
    Defaults to \verb{Inf} which is equivalent to not constrain the re-estimated
    sample size.\cr
    Attention! \verb{max.n} is here \bold{not} a futility criterion like \verb{Nmax}
    in other functions of this package.
  }
  \item{theta0}{
    Assumed ratio of geometric means (T/R) for simulations. If missing, defaults
    to \verb{GMR}.
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
  \emph{et al}.\cr
  The underlying subject data are assumed to be evaluated after log-transformation.
  But instead of simulating subject data, the statistics pe1, mse1 and pe2, SS2 are
  simulated via their associated distributions (normal and
  \ifelse{html}{\out{&chi;<sup>2</sup>}}{\eqn{\chi^{2}}} distributions).
  }
\value{
  Returns an object of class class \verb{"pwrtsd"} with all the input arguments and results
  as components.\cr
  The class class \verb{"pwrtsd"} has a S3 print method.\cr
  The results are in the components:
  \item{pBE}{Fraction of studies found BE.}
  \item{pct_s2}{Percentage of studies continuing to \ifelse{html}{\out{stage&nbsp;2}}{stage~2}.}
  \item{nmean}{Mean of n(total).}
  \item{nrange}{Range (min, max) of n(total).}
  \item{nperc}{Percentiles of the distribution of n(total).}
  \item{ntable}{Object of class \verb{"table"} summarizing the discrete distribution of
  n(total) via its unique values and counts of occurences of these values.\cr
  ntable is only given back if \verb{usePE = FALSE}}
}
\note{
  The computation time is in the magnitude of a few seconds for 100,000 sim\enc{’}{'}s
  on my machine (Intel core i7 2.5 GHz, 12GB RAM) if the non-central \emph{t}
  approximation is used. Thus be a bit patient if you simulate for the Tpye I Error
  'alpha' with 1 Mio sim\enc{’}{'}s.\cr
  Using the crude \verb{pmethod="ls"} on the other hand results in a nearly immediate
  sample size re-estimation.
}
\author{
D. Labes
}
\references{
Golkowski D, Friede T, Kieser M. \emph{Blinded sample size re-estimation in crossover bioequivalence trials.}\cr
Pharm Stat. 2014; 13(3):157--62. \doi{10.1002/pst.1617}

Jones B, Kenward MG. \emph{Design and Analysis of Cross-Over Trials.}\cr
Boca Raton: CRC Press; 3\ifelse{html}{\out{<sup>rd</sup>}}{\eqn{^{rd}}} edition 2014. Chapter 12.

Potvin D, DiLiberti CE, Hauck WW, Parr AF, Schuirmann DJ, Smith RA. \emph{Sequential design approaches for bioequivalence studies with crossover designs.}\cr
Pharm Stat. 2008; 7(4):245--62. \doi{10.1002/pst.294}
}

\seealso{
\code{\link{power.tsd}} for 2-stage studies with interim BE decision.
}
\examples{
# Not run to comply with CRAN policy about examples' run-time;
# minimum number of sim's should be 1E5 for 'power', 1E6 sim's for 'alpha'
\dontrun{
power.tsd.ssr(alpha=0.05, n1=10, GMR=1, CV=0.239, targetpower=0.9,
              pmethod="ls", blind=TRUE, theta0=1.25)
# should give an alpha-inflation 0.072359 (run time <5 seconds)
# repeated with noncentral t-approximation
power.tsd.ssr(alpha=0.05, n1=10, GMR=1, CV=0.239, targetpower=0.9,
              pmethod="nct", blind=TRUE, theta0=1.25)
# should give an alpha-inflation 0.069789 (run time ~20 seconds)
#
# adjusted alpha to control the Type I Error, noncentral t-approx.
power.tsd.ssr(alpha=0.03505, n1=10, GMR=1, CV=0.239, targetpower=0.9,
              pmethod="nct", blind=TRUE, theta0=1.25)
# should control the TIE with 0.049877}
}
