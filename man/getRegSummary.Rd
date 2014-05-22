\name{getRegSummary}
\alias{getRegSummary}
\alias{getRegSummary,GMPM-method}
\alias{getRegSummary,GMPM.mul-method}
\title{Get Regression Summary from GMPM}
\description{
  Creates a summary of the regression analysis from a \code{GMPM} object,
  include parameter estimates and p-values.
}
\usage{
getRegSummary(x)
}
\arguments{
  \item{x}{A \code{GMPM} object}
}
\details{
  No p-values are provided for covariates.
}
\value{

  A list of data.frames.  In the usual case, this will be a list with a
  single element.  For baseline-category multinomial data, there will be
  a single element for each comparison.

}
\seealso{ \code{\link{getMainSummary}} }
