\name{gmpm}
\alias{gmpm}
\alias{gmpmCreate}
\alias{gmpmEstimate}
\alias{gmpmEstimate,GMPM-method}
\title{Creating and Estimating Generalized Multilevel Permutation Models (GMPMs)}
\description{
\code{gmpmCreate} creates a GMPM object without performing permutation runs.
\code{gmpmEstimate} performs permutation runs, and \code{gmpm} is a convenience
function that performs a call to gmpmCreate followed by a call to gmpmEstimate.
}

\usage{

gmpmCreate(formula, family, data=parent.frame(), ivars, gmpmControl=list(), ...)

gmpmEstimate(x, gmpmControl)

gmpm(formula, family, data = parent.frame(), ivars, gmpmControl=list(), ...)

}

\arguments{

  \item{formula}{A symbolic formula describing the regression model to
    be fitted, including a multilevel blocking factor following the '|'
    symbol.  Details of model specification can be found in the section
    \sQuote{Details}.}

  \item{family}{One of the following options: a valid glm class (e.g.,
    binomial, gaussian, poisson, etc.), \dQuote{multinomial}, and
    \dQuote{user} (for a user-defined model).  Options
    \dQuote{multinomial} and \dQuote{user} must be entered within
    quotation marks.}
    
  \item{data}{A data frame containing the data model will be fit to.  If
  left unspecified, gmpm will search for the data within the parent
  environment.}

  \item{ivars}{Names of all independent variables (IVs) included as
  predictors in the model formula, entered as a vector of character
  strings.  All other predictors in the model will be considered
  \sQuote{covariates} and will not be subject to randomization.}

  \item{gmpmControl}{A gmpmControl object used to control the fitting
    function.  This is a list one or more of the following options.  Any
    option not specified in the list will be given a default value.

    \itemize{
      \item{maxruns}{The number of Monte Carlo runs to be performed
      (default = 999).}

      \item{report.interval}{The run interval at which to intermittently
	report, during the fitting process, the number of runs completed and
	the estimated time remaining (default = 10).}

      \item{nCores}{The number of processing cores to use (default =
	"all").  Using more processors will generally result in faster
	performance.  This variable can either be the character
	expression "all" (use all available cores), the character
	expression "all.but.one" (leave one core out so that the system
	is not bogged down) or an integer.  The multicore package must
	be installed, otherwise the number of cores will default to
	one.}
    } }

    \item{x}{A GMPM object}
    
  \item{\dots}{Other arguments to be passed on to the underlying function
    that performs the actual regression (see below).}
}

\details{

  Generalized Multilevel Permutation Models use regression to fit a
  model to the data under the original labeling of experimental
  conditions, and then compares the parameter estimates from this
  \dQuote{original fit} against null-hypothesis distributions obtained
  from fits to random relabelings of the data set.  These functions are
  especially useful for time-series data with difficult or unknown
  sources of dependency.  They can also be used for multivariate
  analyses.  Currently, p-values are provided for categorical predictors
  and interactions between categorical predictors and continuous
  covariates, but not for continuous covariates themselves.

  The \sQuote{formula} option to the model specifies the form of the
  model to be fit.  The left side of the formula identifies the
  dependent variable (DV), which is usually a vector of: real numbers
  (for continuous data), whole numbers (for counts), or 0s and 1s (for
  binomial data).  It can also be a factor with two-levels (binomial
  data) or more (multinomial data).  Lastly, one can also use the
  \code{cbind} syntax for binomial or multinomial DVs (see example
  below).

  The right hand side of the formula specifies the predictors in the
  regression model, using standard R format (see \code{\link{glm}} for
  further information).  The predictors can include both design
  variables (IVs) and covariates.  Any IVs to be randomized must be
  identified in the \sQuote{ivars} option; all other variables are
  assumed to be covariates, and will not be randomized.  All IVs will be
  treated as categorical factors; continuous IVs are not allowed.  The
  IVs are internally coded using \sQuote{deviation} coding, such that
  the codes sum to zero.  The function \code{\link{getFactorCodes}} can
  be used to retrieve this internal coding.

  The \sQuote{family} option specifies the type of dependent variable.
  This option determines which of various regression subroutines are
  used to fit the model.  All of the options for generalized linear
  models (glms) are available (see \code{\link{family}} for details); if
  you choose one of these options then gmpm will call the underlying
  regression function \code{\link{glm}} to perform the actual
  regression.  If the option \sQuote{multinomial} is chosen, then gmpm
  will fit a multinomial model using \code{\link{multinom}} in package
  \pkg{nnet}.  In case tweaks to these underlying functions are
  necessary, it is possible to specify options to \code{gmpm} or
  \code{gmpmCreate} that will be passed along during the fitting process.
  For example, with multinomial models it is often necessary to increase
  the maximum number of iterations, which can be done by the argument
  \sQuote{maxit}.  User-defined models are not yet implemented, but left
  for future development.
}

\value{
  Each function returns a \code{\linkS4class{GMPM}} object.
}

\references{ Barr, D. J. (in preparation).  Generalized Multilevel
  Permutation Models. }

\author{Dale J. Barr <dale.barr@ucr.edu>}

\seealso{

  See example \code{\link{kb07}} for an example of usage with
  binomial data.

  \code{\linkS4class{GMPM}}
}

\examples{
# create a df for a within-subject design
# with no main effect
df1 <- data.frame(SubjID=rep(1:20,each=2),
                  A=factor(rep(c("a1","a2"),times=10)),
                  Y=rnorm(40)*2)

# parametric and non-parametric analyses should look similar
df1.gmpm <- gmpm(Y ~ A | SubjID, gaussian, df1, c("A"))

# parametric and non-parametric analyses should look similar
summary(df1.gmpm)
summary(aov(Y ~ A + Error(SubjID), df1))

# now create horrible dependencies by reduplicating observations
df2 <- rbind(df1,df1,df1,df1,df1,df1,df1,df1,df1)

# parametric analysis more likely to give Type I error
summary(aov(Y ~ A + Error(SubjID), df2))
df2.gmpm <- gmpm(Y ~ A | SubjID, gaussian, df2, c("A"))
}
\keyword{htest}
\keyword{robust}
\keyword{design}
\keyword{ts}
\keyword{nonparametric}
\keyword{multivariate}
