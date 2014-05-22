\name{GMPMSummary-class}
\docType{class}
\alias{GMPMSummary-class}
\alias{initialize,GMPMSummary-method}
\alias{show,GMPMSummary-method}
\title{Class GMPMSummary}
\description{Summary information for a GMPM.}

\section{Objects from the Class}{

  Objects can be created by a call to \code{\link{summary}}.

}

\section{Slots}{
  \describe{

    \item{\code{gmpmInfo}:}{Object of class \code{"list"}; basic
      information about the model.}
    
    \item{\code{gmpmMainSum}:}{Object of class \code{"list"}; the main
      summary of p-values, including multiparameter tests where
      appropriate.}
    
    \item{\code{gmpmRegSum}:}{Object of class \code{"list"}; optional
    summary of the regression, not included by default (see
    \sQuote{showReg} below).}
    
    \item{\code{showReg}:}{Object of class \code{"logical"} whether or
      not the regression summary is to be calculated and displayed (default=TRUE).}
    
  }
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "GMPMSummary")}: ... }
    \item{show}{\code{signature(object = "GMPMSummary")}: ... }
	 }
}
\seealso{
  \code{\link{gmpm}}, \code{\link{GMPM-class}}
}
\keyword{classes}
