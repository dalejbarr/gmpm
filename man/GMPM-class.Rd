\name{GMPM-class}
\docType{class}
\alias{GMPM-class}
\alias{GMPM.glm-class}
\alias{GMPM.mul-class}
\alias{GMPM.user-class}
\alias{initialize,GMPM-method}
\alias{initialize,GMPM.glm-method}
\alias{initialize,GMPM.mul-method}
\alias{initialize,GMPM.user-method}
\alias{show,GMPM-method}
\alias{summary,GMPM-method}

\title{GMPM Model Classes}
\description{S4 objects containing information about GMPMs.}
\section{Objects from the Class}{

  The class \sQuote{GMPM} is a virtual class: No objects may be created
  from it.  The function \code{\link{gmpmCreate}} chooses the type of
  class to be created, either \code{GMPM.glm}, \code{GMPM.mul}, or
  \code{GMPM.user}, depending on the value passed to the \sQuote{family}
  argument of the function \code{\link{gmpmCreate}} or \code{\link{gmpm}}.

  Objects of type \code{GMPM.glm} fit models using the underlying
  regression function \code{\link{glm}}.

  Object of type \code{GMPM.mul} fit models using the underlying function
  \code{\link{multinom}} from package \pkg{nnet}.

  The class \code{GMPM.user} is intended to allow users to specify their
  own regression functions (e.g., calls to nonlinear regression), but
  has not yet been implemented.

}
\section{Slots}{
  \describe{

    \item{\code{df1}:}{a \code{"data.frame"} object containing the model
      frame}
    
    \item{\code{dform}:}{\code{"formula"} specifying the structure of
    the model without the multilevel grouping factor.  This is necessary
    for passing along to the underlying regression function}
    
    \item{\code{mform}:}{The full \code{"formula"} describing the model,
    including the multilevel grouping factor}
    
    \item{\code{munit}:}{Object of type \code{"character"} identifying
    the name of the variable in the model frame that defines the
    multilevel units}
    
    \item{\code{nunits}:}{An integer corresponding to the number of
    multilevel units in the model frame}
    
    \item{\code{gmpmControl}:}{A \code{"list"} with control arguments for
    the fitting function (see \code{\link{gmpm}})}
    
    \item{\code{fitcall}:}{The call to the underlying function, cast as
    a \code{"list"} }
    
    \item{\code{famtype}:}{A human-readable redescription of the family
    of model}
    
    \item{\code{DVname}:}{A string specifying the name of the dependent
    variable}
    
    \item{\code{ivix}:}{A numeric vector giving the indices of the
    permuted variables in the permutation matrix (i.e., variables for
    which p-values are obtained).}
    
    \item{\code{IVinfo}:}{Human-readable information about the
    independent variables (IVs)}
    
    \item{\code{nWithin}:}{Number of \sQuote{within unit} IVs}
    
    \item{\code{nBetween}:}{Number of \sQuote{between unit} IVs}
    
    \item{\code{ivWithin}:}{vector of names of within-unit IVs}
    
    \item{\code{ivBetween}:}{vector of names of between-unit IVs}
    
    \item{\code{ivars}:}{vector of names of all IVs}
    
    \item{\code{IVcoef}:}{A \code{"list"} of all IVs and the names in
    the regression output}
    
    \item{\code{covars}:}{A vector of strings identifying the covariates
    in the model (anything on the right hand side of the formula that is
    not an IV is assumed to be a covariate)}
    
    \item{\code{coefTerms}:}{Object of class \code{"list"} ~~ }
    
    \item{\code{psBetween}:}{A \code{"data.frame"} giving the
    permutation scheme for between-subject variables}
    
    \item{\code{psWithin}:}{A \code{"list"} giving the permutation
    scheme for within-subject variables}
    
    \item{\code{pspace}:}{Number of possible permutations}
    
    \item{\code{pmx}:}{The permutation matrix (or array, for GMPM.mul objects)}
    
    \item{\code{nCellsPerUnit}:}{Number of within-subject cells per
      sampling unit}
    
    \item{\code{ncomp}:}{Number of Monte Carlo runs completed}
  
    \item{\code{convergence}:}{A logical vector specifying whether or
      not a given run converged (\code{GMPM.mul} only)}
    
  }
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "GMPM")}: ... }
    \item{show}{\code{signature(object = "GMPM")}: ... }
    \item{summary}{\code{signature(object = "GMPM")}: ... }
    \item{mdTest}{\code{signature(x = "GMPM", y = "vector")}: ... }
    \item{mdTest}{\code{signature(x = "GMPM", y = "matrix")}: ... }    
	 }
}
\seealso{
  \code{\link{gmpm}}
}
\keyword{classes}
