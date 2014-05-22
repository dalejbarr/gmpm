\name{Mdtest-class}
\docType{class}
\alias{Mdtest-class}
\alias{show,Mdtest-method}
\alias{summary,Mdtest-method}

\title{Class "Mdtest" for Mahalanobis Distance Tests}
\description{An S4 object containing the results of a Mahalanobis
  Distance Test}
\section{Objects from the Class}{
Objects are created by a call to \code{\link{mdTest}}.
}
\section{Slots}{
	 \describe{
    \item{\code{nSections}:}{Number of sections in the results}
    \item{\code{mdsections}:}{list of sections}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "Mdtest")}: ... }
    \item{summary}{\code{signature(object = "Mdtest")}: ... }
	 }
}
\seealso{
  \code{\link{mdTest}}
}
\keyword{classes}
