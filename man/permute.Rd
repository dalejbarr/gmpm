\name{permute}
\alias{permute}
\alias{permute,GMPM,missing-method}
\alias{permute,GMPM,character-method}
\title{Randomly Permute the Model Frame}
\description{
Randomly re-labels the data according to a multilevel relabeling scheme
}
\usage{
permute(x,thisiv)
}
\arguments{
  \item{x}{the \code{GMPM} object to be permuted}
  \item{thisiv}{the independent variable to be tested (or missing)}  
}
\value{
  A permuted version of the object.  If \code{thisiv} is missing, will
  permute all within-unit independent variables.
}
