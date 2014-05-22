\name{getPValue}
\alias{getPValue}
\alias{getPValue,GMPM-method}
\alias{getPValue,matrix-method}
\title{Get p-Values from Permutation Matrix}
\description{
Calculate p-values from the permutation matrix.
}
\usage{
getPValue(x, y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A \code{GMPM} object or a permutation matrix}
  \item{y}{A single number or vector of numbers indexing the columns in
    the permutation matrix for which the calculation is to be performed.}
}
\value{
  A vector of p-values.
}
\seealso{ \code{\link{getNExceeding}} }
