\name{kb07}
\docType{data}
\alias{kb07}
\title{Sample visual world eyetracking data}
\description{
Sample data from a visual-world eyetracking study included to
demonstrate use of gmpm functions on multilevel data with difficult
dependencies.  Data are from a subset of Experiment 2 in Kronmuller
and Barr (2007).  A parametric reanalysis of these data can also be found
in Barr (2008).
}
\usage{
kb07
}
\format{A data frame containing 2464 rows, each one an aggregate over
  multiple trials within the same condition for the same subject.

  \tabular{rlll}{
    [,1] \tab SubjID \tab factor \tab experimental subject\cr
    [,2] \tab Speaker \tab factor \tab identity of the speaker (same,diff)\cr
    [,3] \tab Load \tab factor \tab was listener under cognitive
    load?\cr
    [,4] \tab AggID \tab numeric \tab index of unique
    speaker*load*subject combinations\cr
    [,5] \tab ms \tab factor \tab time in ms from description onset\cr
    [,6] \tab t1 \tab numeric \tab time in 100 ms increments\cr
    [,7] \tab T \tab numeric \tab count of looks to target\cr
    [,8] \tab O \tab numeric \tab count of looks to non-target object\cr
    [,9] \tab X \tab numeric \tab count of frames including looks to
    blank regions or blinks\cr
    [,10] \tab N \tab numeric \tab total number for frames (T+O+X)\cr
    [,11] \tab NT \tab numeric \tab total non-target frames (N-T)
  }
}

\references{ Kronmuller, E. and Barr, D. J. (2007).  Perspective-free
  pragmatics: Broken precedents and the recovery-from-preemption
  hypothesis.  \emph{Journal of Memory and Language}, \bold{56}, 436--455.

  Barr, D. J. (2008).  Analyzing 'visual world' eyetracking data using
  multilevel logistic regression.  \emph{Journal of Memory and Language},
    \bold{59}, 457--474.  }
\author{ Dale J. Barr <dale.barr@ucr.edu> }
\seealso{ \code{\link{gmpm}} }

\examples{
data(kb07)

# first let's fit the data using a binomial model
#
kb07.binom.gmpm <- gmpmCreate(cbind(T,NT) ~ t1*Speaker*Load | SubjID,
                       "binomial", kb07, ivars=c("Speaker","Load"))

# you will need to increase the number of runs in the command below
# (to, say, 999) to get sensible results

kb07.binom.gmpm <- gmpmEstimate(kb07.binom.gmpm, list(maxruns=19))

summary(kb07.binom.gmpm)

# now let's do a more powerful multinomial analysis
# where we break out looks to blank regions and blinks
# into a separate "junk" category (X)

kb07.mnom.gmpm <- gmpmCreate(cbind(O,T,X) ~ t1*Speaker*Load | SubjID,
                     "multinomial", kb07, ivars=c("Speaker","Load"))

# again, you'll need to increase the number of runs to get something
# sensible
kb07.mnom.gmpm <- gmpmEstimate(kb07.mnom.gmpm, list(maxruns=19))

summary(kb07.mnom.gmpm)

}

\keyword{datasets}
