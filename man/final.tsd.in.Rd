\encoding{utf-8}
\name{final.tsd.in}
\alias{final.tsd.in}
\alias{final.2stage.in}
\title{
Analysis after second stage of 2-stage 2x2 crossover design based on
the Inverse Normal method
}
\description{
  Following the design scheme according to \code{power.tsd.in} the function
  performs the analysis after the second stage has been performed.
}
\usage{
final.tsd.in(alpha, weight, max.comb.test = TRUE, GMR1, CV1, n1, df1 = NULL,
             SEM1 = NULL, GMR2, CV2, n2, df2 = NULL, SEM2 = NULL,
             theta1, theta2)
}
\arguments{
  \item{alpha}{
    If one element is given, the overall one-sided significance level (not the
    adjusted level for \ifelse{html}{\out{stage&nbsp;2}}{stage 2}). If two
    elements are given, the adjusted one-sided alpha levels for
    \ifelse{html}{\out{stage&nbsp;1}}{stage 1} and
    \ifelse{html}{\out{stage&nbsp;2}}{stage 2}, respectively.\cr
    If missing, defaults to \verb{0.05}.
  }
  \item{weight}{
    Pre-defined weight(s) of \ifelse{html}{\out{stage&nbsp;1}}{stage 1}.
    Note that using the notation from Maurer et al weight corresponds to
    information fraction, other literature may refer to sqrt(weight) as
    being the weight. \verb{weight} must either contain one element (in case of
    \verb{max.comb.test = FALSE}) or two elements (in case of
    \verb{max.comb.test = TRUE}).\cr
    If missing, defaults to \verb{0.5} for \verb{max.comb.test = FALSE} and to
    \verb{c(0.5, 0.25)} for \verb{max.comb.test = TRUE}.
  }
  \item{max.comb.test}{
    Logical; if \verb{TRUE} (default) the maximum combination test will be
    used, otherwise the standard combination test.
  }
  \item{GMR1}{
    Observed ratio of geometric means (T/R) of
    \ifelse{html}{\out{stage&nbsp;1}}{stage 1} data
    (use \emph{e.g.}, 0.95 for 95\%).
  }
  \item{CV1}{
    Observed coefficient of variation of the intra-subject variability of
    \ifelse{html}{\out{stage&nbsp;1}}{stage 1} (use \emph{e.g.}, 0.3 for 30\%).
  }
  \item{n1}{
    Sample size of \ifelse{html}{\out{stage&nbsp;1}}{stage 1}.
  }
  \item{df1}{
    Optional; Error degrees of freedom of
    \ifelse{html}{\out{stage&nbsp;1}}{stage 1} that can be specified in
    addition to \code{n1}.
  }
  \item{SEM1}{
    Optional; Standard error of the difference of means of
    \ifelse{html}{\out{stage&nbsp;1}}{stage 1} that can be specified in
    addition to \code{CV1}. Must be on additive scale (i.e. usually log-scale).
  }
  \item{GMR2}{
    Observed ratio of geometric means (T/R) of (only)
    \ifelse{html}{\out{stage&nbsp;2}}{stage 2} data
    (use \emph{e.g.}, 0.95 for 95\%).
  }
  \item{CV2}{
    Observed coefficient of variation of the intra-subject variability of (only)
    \ifelse{html}{\out{stage&nbsp;2}}{stage 2} (use \emph{e.g.}, 0.3 for 30\%).
  }
  \item{n2}{
    Sample size of \ifelse{html}{\out{stage&nbsp;2}}{stage 2}.
  }
  \item{df2}{
    Optional; Error degrees of freedom of (only)
    \ifelse{html}{\out{stage&nbsp;2}}{stage 2} that can be specified in
    addition to \code{n2}.
  }
  \item{SEM2}{
    Optional; Standard error of the difference of means of (only)
    \ifelse{html}{\out{stage&nbsp;2}}{stage 2} that can be specified in
    addition to \code{CV2}. Must be on additive scale (i.e. usually log-scale).
  }
  \item{theta1}{
    Lower bioequivalence limit. Defaults to 0.8.
  }
  \item{theta2}{
    Upper bioequivalence limit. Defaults to 1.25.
  }
}
\details{
  The observed values \code{GMR1}, \code{CV1}, \code{n1} must be obtained
  using data from stage 1 only, and \code{GMR2}, \code{CV2}, \code{n2} must
  be obtained using data from stage 2 only. This may be done via the usual
  ANOVA approach.\cr\cr
  The optional arguments \code{df1}, \code{SEM1}, \code{df2} and \code{SEM2}
  require a somewhat advanced knowledge (provided in the raw output from for
  example the software \acronym{SAS}, or may be obtained via \code{emmeans::emmeans}).
  However, it has the advantage that if there were missing data the exact
  degrees of freedom and standard error of the difference can be used,
  the former possibly being non-integer valued (e.g. if the
  Kenward-Roger method was used).
}
\value{
  Returns an object of class \verb{"evaltsd"} with all the input arguments and results
  as components. As part of the input arguments a component \code{cval} is also
  presented, containing the critical values for \ifelse{html}{\out{stage&nbsp;1}}{stage 1} and 2 according to the
  input based on \code{alpha}, \code{weight} and \code{max.comb.test}.\cr
  The class \verb{"evaltsd"} has an S3 print method.\cr\cr
  The results are in the components:
  \item{z1}{Combination test statistic for first null hypothesis (standard
  combination test statistic in case of \code{max.comb.test = FALSE} or
  maximum combination test statistic in case of \code{max.comb.test = TRUE})}
  \item{z2}{Combination test statistic for second null hypothesis (standard
  combination test statistic in case of \code{max.comb.test = FALSE} or
  maximum combination test statistic in case of \code{max.comb.test = TRUE})}
  \item{RCI}{Repeated confidence interval for
  \ifelse{html}{\out{stage&nbsp;2}}{stage 2}.}
  \item{MEUE}{Median unbiased point estimate as estimate for the final adjusted
  geometric mean ratio after stage 2.}
  \item{stop_BE}{Logical, indicating whether \acronym{BE} can be concluded after
  \ifelse{html}{\out{stage&nbsp;2}}{stage 2} or not.}
}
\author{
B. Lang
}
\references{
König F, Wolfsegger M, Jaki T, \enc{Schütz}{Schuetz} H, Wassmer G.\cr
\emph{Adaptive two-stage bioequivalence trials with early stopping and sample size re-estimation.}\cr
Vienna: 2014; 35\ifelse{html}{\out{<sup>th</sup>}}{\eqn{^{th}}} Annual Conference of the International Society for Clinical Biostatistics. Poster P1.2.88 \cr
\doi{10.13140/RG.2.1.5190.0967}.

Patterson SD, Jones B. \emph{Bioequivalence and Statistics in Clinical Pharmacology.}\cr
Boca Raton: CRC Press; 2\ifelse{html}{\out{<sup>nd</sup>}}{\eqn{^{nd}}} edition 2017.

Maurer W, Jones B, Chen Y. \emph{Controlling the type 1 error rate in two-stage
sequential designs when testing for average bioequivalence.}\cr
Stat Med. 2018; 37(10): 1587--1607. \doi{10.1002/sim.7614}.

Wassmer G, Brannath W. \emph{Group Sequential and Confirmatory Adaptive Designs
in Clinical Trials.}\cr
Springer 2016. \doi{10.1007/978-3-319-32562-0}.
}
\seealso{
\code{\link{power.tsd.in}, \link{interim.tsd.in}}
}
\examples{
# Example from Maurer et al.
final.tsd.in(GMR1 = exp(0.0424), CV1 = 0.3682, n1 = 20,
             GMR2 = exp(-0.0134), CV2 = 0.3644, n2 = 36)
# Example 2 from Potvin et al.
final.tsd.in(GMR1 = 1.0876, CV1 = 0.18213, n1 = 12,
             GMR2 = 0.9141, CV2 = 0.25618, n2 = 8)
}
