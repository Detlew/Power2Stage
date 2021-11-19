\encoding{utf-8}
\name{sampleN2.TOST}
\alias{sampleN2.TOST}
\title{
Sample size re-estimation of adaptive 2-stage BE studies in 2x2 crossover and 
parallel designs based on power of TOST
}
\description{
  This function estimates the necessary sample size of
  \ifelse{html}{\out{stage&nbsp;2}}{stage 2} to have at least a given power.
  }
\usage{
sampleN2.TOST(alpha = 0.0294, CV, n1, theta0 = 0.95,
              theta1 = 0.8, theta2 = 1.25, targetpower = 0.8,
              design = "2x2", method = "exact", imax = 100)
}
\arguments{
  \item{alpha}{
    Alpha value for the final analysis of pooled data. Defaults to Pocock\enc{’}{'}s
    alpha setting \verb{alpha=0.0294}.
  }
  \item{CV}{
    Coefficient of variation of the intra-subject variability as ratio.
  }
  \item{n1}{
    Sample size of \ifelse{html}{\out{stage&nbsp;1}}{stage 1}.
  }
  \item{theta0}{
    True ratio of T/R for simulating. Defaults to 0.95 argument if missing.
  }
  \item{theta1}{
    Lower bioequivalence limit. Defaults to 0.8.
  }
  \item{theta2}{
    Upper bioequivalence limit. Defaults to 1.25.
  }
  \item{targetpower}{
    Power to achieve at least. Must be >0 and <1.
  }
  \item{design}{
    Character string describing the study design. Implemented are \verb{"2x2"}
    and \verb{"parallel"}.
  }
  \item{method}{
    Method for calculation of the power. Implemented are \verb{"exact"}
    (exact calculation via Owen\enc{’}{'}s Q), \verb{"nct"} (approximate
    calculation via non-central \emph{t}-distribution, and \verb{"shifted"}
    (approximate calculation via shifted central \emph{t}-distribution
    like in the paper of Potvin \emph{et al.}\cr
    Defaults to \verb{"exact"}.
  }
  \item{imax}{
    Maximum number of steps in sample size search.\cr
    Defaults to 100. Adaption only in rare cases needed.
  }
}
\details{
  The sample size is calculated via iterative evaluation of power of the
  TOST procedure.\cr
  Start value for the sample size search is taken from a large sample
  approximation according to Zhang, modified.
  }
\value{
  A data.frame with the input and results will be returned.\cr
  The "Sample size" column contains the sample size for the second stage.
  }
\note{
  Since in the final analysis one degree of freedom less than in a fixed sample
  design (or in \ifelse{html}{\out{stage&nbsp;1}}{stage 1}) is used, power will
  be slightly lower than the one obtained with function \verb{sampleN.TOST} of
  package \verb{PowerTOST}. However, different \emph{sample sizes} are extremely
  unlikely.
}
\references{
Potvin D, DiLiberti CE, Hauck WW, Parr AF, Schuirmann DJ, Smith RA. \cr
\emph{Sequential design approaches for bioequivalence studies with crossover designs.}\cr
Pharm Stat. 2008; 7(4):245--62. \doi{10.1002/pst.294}

Zhang P. \emph{A Simple Formula for Sample Size Calculation in Equivalence Studies.}\cr
J Biopharm Stat. 2003; 13(3):529--38. \doi{10.1081/BIP-120022772}
}
\author{
H. \enc{Schütz}{Schuetz}, D. Labes
}
\examples{
# using all the defaults, CV of 25\% and 12 subjects in stage 1
print(sampleN2.TOST(CV=0.25, n1=12), row.names=FALSE)
# should give a stage 2 sample size of 22 and achieved power ~0.812
# CV 10\% and 12 subjects in stage 1
print(sampleN2.TOST(CV=0.1, n1=12), row.names=FALSE)
# should give a sample size of 0 (second stage not reasonable
# since power ~0.973 was already achieved in stage 1)
}
