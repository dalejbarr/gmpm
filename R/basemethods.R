setClass("GMPM",
         representation(
                        df1="data.frame", # the model.frame
                        dform="formula", # design formula (including covars)
                        mform="formula", # full model formula
                        munit="character", # multilevel sampling unit
                        nunits="numeric", # number of units sampled
                        gmpmControl="list", # control of fitting functions
                        fitcall="list", # call to fitting function
                        famtype="character", # type of data
                        DVname="character", # name of DV
                        ivix="numeric", # indices of IVs
                        IVinfo="data.frame", # info about IVs
                        nWithin="numeric", # nWithin unit variables
                        nBetween="numeric", # nBetween unit variables
                        ivWithin="character", # list of within IVs
                        ivBetween="character", # list of between IVs
                        minN="numeric", # nObs in smallest cell in design
                        ivars="vector", # list of names of IVs
                        IVcoef="list", # names of factor vars in glm output
                        covars="character", # list of names of covars
                        coefTerms="list", # names of variables from fit output
                                          # w/interactions separated.

                        psBetween="list", # permn scheme (Betw unit IVs)
                        psWithin="list", # permutation scheme (Within unit IVs)
                        nwrep="numeric", # n. within reps per withinIV
                        pspace="numeric", # size of permutation space
                        nSections="numeric", # number of permutation sections for estimation
                        psec="list", # permutation sections
                        
                        pmx="matrix", # matrix of permutation coefficients
                        nCellsPerUnit="numeric", # nCells per sampling unit

                        ncomp="numeric", # n of runs completed
                        ndigits="numeric", # n of digits to round p value to
                        nCores="integer", # N of available processing cores
                        "VIRTUAL"), # factor matrix for the model
         
         prototype=prototype(
           nunits=1, nWithin=0, nBetween=0, ncomp=0),
         )

setClass(Class="GMPMSummary",
         representation(
                        gmpmInfo="list", # misc. info about gmpm object
                        gmpmMainSum="list", # list of data frames
                                           # with summary info
                        gmpmRegSum="list", # main regression
                        showReg="logical" # whether to show reg coef?
                        ),
         prototype(showReg=TRUE)
         )

setClass(
         Class="GMPM.glm",
         representation(
                        coef0="numeric", # vector of original coefficients
                        family="list"
                        ),
         contains="GMPM"
         )

setClass(
         Class="GMPM.mul",
         representation(
                        coef0="matrix", # vector of original coefficients
                        family="character",
                        pmx="array",
                        convergence="vector" # did it converge?
                        ),
         prototype(family="multinomial",famtype="multinomial"),
         contains="GMPM"
         )

setClass(
         Class="GMPM.user",
         representation(
                        family="character"
                        ),
         prototype(family="user",famtype="user"),
         contains="GMPM"
         )

setMethod("initialize",
          signature(.Object = "GMPM"),
          function (.Object,
                    formula, family, data,
                    ivars, gmpmControl)
          {
#            print(">>>> initializing (GMPM)")            
            return(.Object)
          }
          )

setMethod("initialize",
          signature(.Object="GMPM.glm"),
          function(.Object, family=gaussian, ...)
          {
#            print(">>>> initializing (GMPM.glm)")
            if (is.character(family)) {
              family <- get(family, mode = "function",
                            envir = globalenv())
            } else {}

            if (is.function(family)) {
              family <- family()
            } else {}

            if (is.null(family)) {
              print(family)
              stop("'family' not recognized")
            } else {}
            
            ff <- family
            if (class(family) == "family") {
              fname <- ff[[1]]
              ltype <- ff[[2]]
            } else {
              fname <- deparse(substitute(family))
              ltype <- as.list(ff)[[1]]
            }

            .Object@famtype <- paste(fname,
                                     "(link=", ltype,
                                     ")", sep="")

#            print(.Object@famtype)

            .Object@family <- unclass(family)

            return(.Object)
          }
          )

setMethod("initialize",
          signature(.Object="GMPM.mul"),
          function(.Object, ...)
          {
#            print(">>>> initializing (GMPM.mul)")
#            callNextMethod()
            require(nnet)
            return(.Object)            
          }
          )

setMethod("initialize",
          signature(.Object="GMPM.user"),
          function(.Object, ...)
          {
#            print(">>>> initializing (GMPM.user)")
            cat("Warning: User must supply fitting function (see ?createCall for details).\n")
            return(.Object)
          }
          )

setMethod("initialize",
          signature(.Object="GMPMSummary"),
          function(.Object, gmpmInfo, gmpmMainSum=NULL, gmpmRegSum=NULL)
          {
#            print(">>>> initializing (GMPMSummary)")
            .Object@gmpmInfo <- gmpmInfo
            if (!is.null(gmpmMainSum)) {
              .Object@gmpmMainSum <- gmpmMainSum
            } else {}
            if (!is.null(gmpmRegSum)) {
              .Object@gmpmRegSum <- gmpmRegSum
            } else {}
            return(.Object)
          }
          )

setMethod("show",
    signature(object = "GMPM"),
    function (object) 
    {
      xsum <- summary(object)
      print(xsum)
      return()
    }
)

#setMethod("coef",
#    signature(object = "GMPM"),
#    function (object) 
#    {
#      return(gmpCoef(object))
#    }
#)

#setMethod("coefficients",
#    signature(object = "GMPM"),
#    function (object) 
#    {
#      return(gmpCoef(object))
#    }
#)

setMethod("show",
          signature(object = "GMPMSummary"),
          function(object)
          {
            cat("\n")
            x <- object@gmpmInfo
            
            if (x$nunits == 1) {
              cat("Single-level data with", x$nobs, "observations.\n\n")
            } else {
              cat("Multilevel data with ", x$nobs,
                  "observations \nover ", x$nunits,
                  "units,", paste("defined by '", x$munit, "'.", sep=""),
                  "\n\n")
            }
            cat("Dependent variable", paste("'", x$DVname, "'", sep=""),
                "of type", x$famtype, "\n")

            # something for multinomial data here.
            if (x$famtype == "multinomial") {
              cat("Levels: ", x$DVlevels[1], " (baseline) ",
                  "versus ", x$DVlevels[2:length(x$DVlevels)],"\n")
            }
            cat("\n")
            
            cat("Independent variables:\n")
            print(x$IVinfo)
            cat("\n")

            if (length(x$covars) > 0) {
              cat("Covariates: ", paste(x$covars, collapse=", "), "\n\n")
            }
            
            cat("Model: ")
            print(x$mform)
            cat("\n")

            if (!object@showReg) {
              if (is.matrix(x$coef0)) {
                dft <- as.data.frame(cbind(rownames(t(x$coef0)),
                                           round(t(x$coef0),4)))
                rownames(dft) <- 1:length(rownames(dft))
                colnames(dft) <- c("Coefficient",
                                   paste(rownames(x$coef0),
                                         "vs",x$DVlevels[1], sep=" "))
              } else {
                dft <- data.frame(Coef=names(x$coef0), Estimate=round(x$coef0,4))
                rownames(dft) <- 1:length(x$coef0)
                colnames(dft) <- c("Coefficient", "Estimate")
              }
              print(dft)
            } else {
              x <- object@gmpmRegSum
              if (length(x) > 0) {
                cat("Summary of Individual Regression Parameters:\n")
                if (length(x) == 1) {
                  print(x[[1]])
                } else {
                  snames <- names(x)
                  for (q in 1:length(x)) {
                    if (q > 1) {
                      cat("\n")
                    } else {}
                    cat("------ ", snames[q], " ------\n")
                    print(x[[q]])
                  }
                }
                cat("---\n")
                cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
                cat("\n")
              } else {}
            }

            cat("\n")

                                        # now come the main results
            mainSum <- object@gmpmMainSum
            nSections <- length(mainSum)
            if (nSections > 0) {
              cat(">>>>>>>>> SUMMARY OF MAIN RESULTS <<<<<<<<<\n\n")
              secnames <- names(mainSum)
              for (i in 1:nSections) {
                x <- mainSum[[i]]
                if (dim(x)[1] > 0) {
                  if (i > 1) {
                    cat("\n")
                  } else {}
                  if (nSections > 1) {
                    cat("-----", secnames[i], "-----\n")
                  }
                  print(x)
                } else {}
              }
              cat("---\n")
              cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")

            }
            cat("\n")

            if (object@gmpmInfo$ncomp > 1) {
              cat("All p-values based on", object@gmpmInfo$ncomp,
                  "Monte Carlo samples\n\n")
              #"from ", object@gmpmInfo$pspace,
              #    "possible permutations.\n\n")
              
              if (object@gmpmInfo$ncomp < 999) {
                cat("Warning: Too few Monte Carlo samples for reliable p-values.\n", "Consider increasing 'maxruns'.\n")
              } else {}
            }
          }
          )

setMethod("summary",
    signature(object = "GMPM"),
          function (object, showReg=TRUE, ...) 
          {
#            print("~~~ in summary (GMPM) ~~~")
            x <- object

            gmpmInfo <- list()
            gmpmInfo$nunits <- x@nunits
            gmpmInfo$nobs <- dim(x@df1)[1]
            gmpmInfo$munit <- x@munit
            gmpmInfo$DVname <- x@DVname
            gmpmInfo$famtype <- x@famtype
            gmpmInfo$IVinfo <- x@IVinfo
            gmpmInfo$mform <- x@mform
            gmpmInfo$ncomp <- x@ncomp
            gmpmInfo$pspace <- x@pspace
            gmpmInfo$coef0 <- x@coef0
            gmpmInfo$covars <- x@covars
            if (x@famtype == "multinomial") {
              gmpmInfo$DVlevels <- .getDVlevels(x)
            }

            if (x@ncomp <= 1) {
              xsum <- new("GMPMSummary",
                          gmpmInfo)
              xsum@showReg <- FALSE
            } else {
                                        # build main summary.


              ###########################
              # loop for each row in permutation array?
              # in order to report results for multinomial
              ###########################
              
                                        # build regression summary.
              if (showReg) {
                gmpmRegSum <- getRegSummary(object)
              } else {
                gmpmRegSum <- data.frame()
              }

              gmpmMainSum <- list()
              faclist <-
                attr(attr(x@df1,"terms"),"factors")[-1,]
              if (is.vector(faclist)) {
                allvars <- x@ivars
                nTests <- 1
              } else {
                allvars <- rownames(faclist)
                nTests <- dim(faclist)[2]
              }
              if (length(x@covars) > 0) {
                m <- match(x@covars, colnames(faclist))
                m <- m[!is.na(m)]
                faclist <- faclist[,-m]
              }
              if (is.vector(faclist)) {
                nTests <- 1
              } else {
                nTests <- dim(faclist)[2]
              }
              # build main summary
              gmpmMainSum <-
                getMainSummary(x)

              xsum <- new("GMPMSummary",
                          gmpmInfo, gmpmMainSum, gmpmRegSum)
              xsum@showReg = showReg

            }
            
#            print("... exiting summary (GMPM) ...")
            return(xsum)
          }
          )
