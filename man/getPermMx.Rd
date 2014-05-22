\name{getPermMx}
\alias{getPermMx}
\alias{getPermMx,GMPM-method}
\title{Get the Permutation Matrix from a GMPM object}
\description{
Retrieves the permutation matrix from a GMPM object.
}
\usage{
getPermMx(x)
}
\arguments{
  \item{x}{A \code{GMPM} object}
}

\value{
  A matrix object.  The first row of the permutation matrix includes the parameter
  estimates from the original fit.  Each remaining row has estimates
  from individual Monte Carlo runs.
}
