###### first come the helper functions

.getMinN <- function(df1,ivs) {
  # find the cell with the minimum no. of observations
  f1 <- combFreqs(df1,ivs)
  return(min(f1$N))
}

sortDF <- function(x, ovec, fixRowNames=TRUE) {
  slst <- list()
  for (i in 1:length(ovec)) {
    slst[[i]] <- x[,ovec[i]]
  }

  df1 <- x[do.call(order, slst),]

  if (fixRowNames) {
    rownames(df1) <- 1:nrow(df1)
  } else {}
  
  return(df1)
}

combFreqs <- function(x, ovec) {

  if (length(ovec)==0) {
    x$N <- 1
    return(x)
  }
  
  N <- 1:nrow(x)
  bylist <- list()
  for (i in 1:length(ovec)) {
    bylist[[i]] <- x[,ovec[i]]
  }
  names(bylist) <- ovec
  
  xt <- aggregate(list(N=N), by=bylist, length)
  return(sortDF(xt,ovec))
}

getDFix <- function(x, keys) {
  # return row indices for rows matching variable-value pairs in keys
  ex1 <- rep(NA, length(keys))
  varnames <- names(keys)
  for (i in 1:length(keys)) {
    if (is.factor(keys[[i]]) || is.character(keys[[i]])) {
      str1 <- "'"
    } else {
      str1 <- ""
    }
    ex1[i] <- paste("(x[,'",varnames[i],"']==",str1,keys[[i]],str1,")",sep="")
  }
  e1 <- parse(text=paste(ex1, collapse=" & "))
  return((1:nrow(x))[eval(e1,list(x=x))])
}

permuteNA <- function(x) {
  s1 <- sample(x)
  return(c(s1[!is.na(s1)],s1[is.na(s1)]))
}

###### now come the S4 methods

setMethod(".initFinal",
          signature(object="GMPM"),
          function(object) {            
            return(object)
          })

setMethod(".initFinal",
          signature(object="GMPM.mul"),
          function(object) {

          # make sure that if it is a single variable, it is coded as a factor
            if (length(grep("cbind", object@DVname))==0) {
              if (!is.matrix(object@df1[,object@DVname])) {
                if (!is.factor(object@df1[,object@DVname])) {
                  object@df1[,object@DVname] <- factor(object@df1[,object@DVname])
                  warning("Converting '", object@DVname, "' to a factor")
                } else {}
              } else {}
            } else {}
            
            return(object)
          })

setMethod("getModelFrame",
          signature(object="GMPM"),
          function(object) {
            return(object@df1)
          })

setMethod("getFactorCodes",
          signature(object="GMPM"),
          function(object) {
            ivs <- object@ivars
            nIVs <- length(object@ivars)
            fcodes <- list()
            for (i in 1:nIVs) {
              mx <- attr(object@df1[,ivs[i]],"contrasts")
              colnames(mx) <- object@IVcoef[[ivs[i]]]
              fcodes[[ivs[i]]] <- mx
            }
            return(fcodes)
          })

setMethod(".getFactorLabelsFromFit",
          signature(object="GMPM"),
          function(object, ivar) {
            fcall <- as.list(object@fitcall)
            dform <- object@dform
            lhs <- strsplit(deparse(dform), "~")[[1]][1]
            for (i in 1:length(object@ivars)) {
              nform <- paste(lhs, "~", object@ivars[i], sep="")
              fcall$formula <- nform
              if (object@famtype == "multinomial") {
                capture.output(nn1 <- colnames(coef(eval(as.call(fcall)))))
              } else {
                nn1 <- names(coef(eval(as.call(fcall))))
              }
              object@IVcoef[[object@ivars[i]]] <- nn1[2:length(nn1)]
            }
            return(object@IVcoef)
          })

setMethod(".getPredictorsFromFaclist",
          signature(x="GMPM"),
          function(x, faclist, allvars, nTests, j) {
            if (nTests > 1) {
              ivinc <- allvars[faclist[,j]==1]
            } else {
              ivinc <- allvars
            }

            setBetw <- intersect(ivinc,x@ivBetween)
            if (length(setBetw) > 0) {
              pmx <- (1:length(x@ivBetween))[x@ivBetween==setBetw[1]]
            } else {
              pmx <- length(x@psec)
            }
            
            if (length(intersect(ivinc, x@covars)) > 1) {
                                        # there is more than one co-variate in the term.
                                        # we need to perform a union rather than an intersection
              cvartmp <- intersect(ivinc, x@covars)
              cvars <- c()
              ivartmp <- intersect(ivinc, x@ivars)
              for (i in 1:length(cvartmp)) {
                cvars <- c(cvars, .getIXfromIV(x, c(ivartmp, cvartmp[i]), TRUE))
              }
            } else {
              cvars <- .getIXfromIV(x, ivinc, TRUE)
            }

            return(list(pmx=pmx, cvars=cvars))
          })

setMethod("testDV",
          signature(x="GMPM.mul"),
          function(x, excludeLevels, byCovar=FALSE) {
            DVlevels <- .getDVlevels(x)
            if (!missing(excludeLevels)) {
              if (is.character(excludeLevels)) {
                excludeLevels <- c(excludeLevels)
              } else {}
              if (DVlevels[1] %in% excludeLevels) {
                stop("Can't exclude baseline level '", DVlevels[1], "' from analysis.")
              } else {}
              ltest <- setdiff(DVlevels[2:length(DVlevels)], excludeLevels)
              if (length(ltest) < 2) {
                stop("Only ", length(ltest), " regions to test.\nMinimum of 2 required.")
              }
            } else {
              ltest <- DVlevels[2:length(DVlevels)]
            }
            ff <- .prepareMainSum(x, byCovar)
            nTests <- ff$nTests
            nDiffs <- length(ltest)-1
            pwid <- dim(x@pmx)[3]
            mx <- matrix(nrow=dim(x@pmx)[1], ncol=nDiffs*pwid)
            colnames(mx) <- rep(dimnames(x@pmx)$coef, nDiffs)
            newDVname <- paste("(",
                               paste(c(DVlevels[1], ltest), collapse=",", sep=""),
                               ")", sep="")
            for (k in 1:nDiffs) {
              ix0 <- (k-1)*(pwid)+1
              ix1 <- ix0+pwid-1
              mx[,ix0:ix1] <- x@pmx[,ltest[1],]-x@pmx[,ltest[k+1],]
            }
            mlist <- list()
            for (j in 1:nTests) {
              cvars <- .getPredictorsFromFaclist(x, ff$faclist, ff$allvars, ff$nTests, j)
              ctest <- cvars + rep(rep(0:(nDiffs-1))*pwid, each=length(cvars))
              tname <- paste(colnames(ff$faclist)[j], ":", newDVname, sep="")
              mlist[[tname]] <- ctest
            }
            return(mdTest(mx, mlist))
          })

setMethod(".getDVlevels",
          signature(x="GMPM.mul"),
          function(x) {
            if(length(grep("cbind", x@DVname))>0) {
              f1 <- strsplit(x@DVname, "\\(")[[1]][2]
              f2 <- gsub(")", "", f1)
              f3 <- gsub(" ", "", f2)
              return(strsplit(f3, ",")[[1]])
            }
            
            if (is.matrix(x@df1[,x@DVname]))
              return(colnames(x@df1[,x@DVname]))
            else if (is.factor(x@df1[,x@DVname]))
              return(levels(x@df1[,x@DVname]))
            else {}

            stop("unsure of how DV '",x@DVname,"' is represented\n",
                 "(was not factor, matrix, or cbind.)")
          })

setMethod(".writeFit",
          signature(x="GMPM"),
          function(x, y, outfile, append) {
            cat(y, "\n", file=outfile, append=append)
          })

setMethod(".writeFit",
          signature(x="GMPM.mul"),
          function(x, y, outfile, append) {
            nRows <- dim(y)[1]
            for (i in 1:nRows) {
              cat(y[i,], " ", file=outfile, append=append)
            }
            cat("\n", file=outfile, append=TRUE)
          })

setMethod("appendToPmx",
          signature(x="GMPM", y="GMPM"),
          function(x, y) {
            nameObject <- deparse(substitute(x))            
            if (length(dim(x)) != length(dim(y))) {
              cat("GMPM source object permutation matrix has dimensions ", dim(y),
                  "\n")              
              cat("GMPM destination object permutation matrix has dimensions ", dim(x),
                  "\n")
              stop("These GMPM objects do not look the same.")
            }
            pmxSrc <- getPermMx(y)[-1,]
            pmxThis <- getPermMx(x)
            if (dim(pmxThis)[1]==0) {
              x@pmx <- getPermMx(y)
            } else {
              x@pmx <- rbind(pmxThis, pmxSrc)
            }
            cat("appended ", length(pmxSrc[,1]), " rows\n")            
            x@ncomp <- dim(x@pmx)[1]-1
            warning("Error checking not implemented yet; \nPlease ensure these two GMPM objects have the same underlying model / data.")
            assign(nameObject, x, envir=parent.frame())            
            return(invisible(x))
          })

setMethod("appendToPmx",
          signature(x="GMPM", y="character"),
          function(x, y) {
            nameObject <- deparse(substitute(x))            
            ff <- read.table(y)
            pmxSrc <- as.matrix(ff[-1,])
            rownames(pmxSrc) <- NULL
            pmxThis <- getPermMx(x)
            if (dim(pmxThis)[1]==0) {
              pmxSrc <- as.matrix(ff)
              colnames(pmxSrc) <- colnames(x@coef0)
              x@pmx <- pmxSrc
            } else {
              colnames(pmxSrc) <- colnames(pmxThis)
              x@pmx <- rbind(pmxThis, pmxSrc)
            }
            x@ncomp <- dim(x@pmx)[1]-1
            cat("appended ", length(pmxSrc[,1]), " rows\n")
            
            assign(nameObject, x, envir=parent.frame())            
            return(invisible(x))
          })

setMethod("appendToPmx",
          signature(x="GMPM.mul", y="GMPM"),
          function(x, y) {
            nameObject <- deparse(substitute(x))            
            if (length(dim(x)) != length(dim(y))) {
              cat("GMPM source object permutation matrix has dimensions ", dim(y),
                  "\n")              
              cat("GMPM destination object permutation matrix has dimensions ", dim(x),
                  "\n")
              stop("These GMPM objects do not look the same.")
            }
            pmxSrc <- getPermMx(y)[-1,,]
            srclen <- dim(pmxSrc)[1]
            pmxThis <- getPermMx(x)
            destlen <- dim(pmxThis)[1]
            nDVlevels <- dim(x@coef0)[1]
            nCoef <- dim(x@coef0)[2]
            dmn <- c("run","dv","coef")
            if (destlen > 0) {
              x@pmx <- array(dim=c(srclen+destlen,nDVlevels,nCoef),
                             dimnames=dmn)
              x@pmx[1:destlen,,] <- pmxThis
              x@pmx[(destlen+1):(srclen+destlen),,] <- pmxSrc
            } else {
              pmxSrc <- getPermMx(y)
              srclen <- dim(pmxSrc)[1]
              x@pmx <- getPermMx(y)
            }
            cat("appended ", srclen, " rows\n")            
            x@ncomp <- dim(x@pmx)[1]-1
            warning("Error checking not implemented yet; \nPlease ensure these two GMPM objects have the same underlying model / data.")
            assign(nameObject, x, envir=parent.frame())            
            return(invisible(x))
          })

setMethod("appendToPmx",
          signature(x="GMPM.mul", y="character"),
          function(x, y) {
            nameObject <- deparse(substitute(x))            
            ff <- read.table(y)
            pmxSrc <- as.matrix(ff[-1,])
            rownames(pmxSrc) <- NULL
            pmxThis <- getPermMx(x)
            colnames(pmxSrc) <- colnames(pmxThis)
            x@pmx <- rbind(pmxThis, pmxSrc)
            x@ncomp <- dim(x@pmx)[1]-1
            cat("appended ", length(pmxSrc[,1]), " rows\n")
            
            assign(nameObject, x, envir=parent.frame())            
            return(invisible(x))
          })

setMethod(".prepareMainSum",
          signature(x="GMPM"),
          function(x, byCovar=FALSE) {

            # figure out tests we need to run from model frame
            faclist <-
              attr(attr(x@df1,"terms"),"factors")[-1,]
            if (is.vector(faclist)) {
              allvars <- x@ivars
              nTests <- 1
            } else {
              allvars <- rownames(faclist)
              nTests <- dim(faclist)[2]
            }

            # do not test main effects of covars
            if (length(x@covars) > 0) {
              m <- match(x@covars, colnames(faclist))
              m <- m[!is.na(m)]
              faclist <- faclist[,-m]
            }

            if ((!byCovar) && (length(x@covars) > 1)) {
              fl1 <- faclist[x@covars[-1],]
              if (is.vector(fl1)) {
                test1 <- fl1==0    
              } else {
                test1 <- colSums(fl1)==0    
              }
              newfl <- faclist[,faclist[x@covars[1],]==1 | test1]
              r1 <- newfl[x@covars[1],]
              c1 <- names(r1)[r1==1]
              newfl[x@covars[-1],c1] <- 1
              flnames <- strsplit(colnames(newfl), ":")
              srep <- paste("(",paste(x@covars, collapse=","),")",sep="")
              ncn <- rep("", length(flnames))
              for (i in 1:length(flnames)) {
                flnames[[i]][flnames[[i]]==x@covars[1]] <- srep
                ncn[i] <- paste(flnames[[i]],collapse=":",sep="")
              }
              colnames(newfl) <- ncn
              faclist <- newfl
            }
            
            if (is.vector(faclist)) {
              nTests <- 1
            } else {
              nTests <- dim(faclist)[2]
            }
            
            return(list(faclist=faclist,
                        allvars=allvars,
                        nTests=nTests))
          })

setMethod(".regSumProc",
          signature(x="GMPM"),
          function(x, psec, index) {
            #print("~~~ in .regSumProc ~~~")
            coef0 <- psec[[1]][1,]

            gg <- .prepareMainSum(x,FALSE)
            psecmap <- rep(NA, length(coef0))
            names(psecmap) <- names(coef0)
            for (i in 1:gg$nTests) {
              hh <- .getPredictorsFromFaclist(x,gg$faclist,gg$allvars,gg$nTests,i)
              for (j in 1:length(hh$cvars)) {
                psecmap[hh$cvars[j]] <- hh$pmx
              }
            }

            c2 <- coef0
            se <- rep(NA, length(c2))
            nexceed <- rep(NA, length(c2))
            pval <- rep(NA, length(c2))              
            for (i in 1:length(c2)) {
              if (!is.na(psecmap[i])) {
                pmx <- psec[[psecmap[i]]]
                se[i] <- sd(pmx[,i])
                nexceed[i] <- getNExceeding(pmx, i)
                pval[i] <- getPValue(pmx, i)
              } else {}
            }

            c2 <- round(c2,4)
            se <- round(se,4)
            pval <- round(pval,x@ndigits)
            gmpmRegSum <- data.frame(Coef=names(c2),
                                    Estimate=c2, se, nexceed,
                                    pval, .getSig(pval))
            rownames(gmpmRegSum) <- 1:length(c2)
            colnames(gmpmRegSum) <- c("Coefficient", "Estimate",
                                     "Std. Error", "N>=orig", "p-value", " ")
            return(gmpmRegSum)
          })

setMethod("getRegSummary",
          signature(x="GMPM"),
          function(x) {
            ff <- list(.regSumProc(x, x@psec, x@ivix))
            names(ff) <- "Main Regression"
            return(ff)
          })

setMethod("getRegSummary",
          signature(x="GMPM.mul"),
          function(x) {
            #print("~~~ in getRegSummary (GMPM.mul) ~~~")
            mlist <- list()
            DVlevels <- .getDVlevels(x)
            nDVlevels <- dim(x@coef0)[1]            
            mnames <- paste(DVlevels[2:length(DVlevels)], DVlevels[1],
                  sep=" versus ")
            for (i in 1:nDVlevels) {
              psec <- .collapseMultinomPmx(x, i)
              mlist[[i]] <- .regSumProc(x, psec, x@ivix)
            }
            names(mlist) <- mnames
            #print("... exiting getRegSummary (GMPM.mul) ...")
            return(mlist)
            
          })

setMethod(".mainSumProc",
          signature(x="GMPM"),
          function(x, faclist, allvars, nTests, psec) {
            #print("~~~ in .mainSumProc ~~~")

            #coef0 <- pmx[1,]
            nge <- rep(NA, nTests)
            pval <- rep(NA, nTests)
            mcoef <- rep(NA, nTests)

            for (j in 1:nTests) {
              hh <- .getPredictorsFromFaclist(x, faclist, allvars, nTests, j)
              pmx <- psec[[hh$pmx]]
              cvars <- hh$cvars
              
              if (length(cvars) > 1) {
                mdt <- mdTest(pmx, cvars)
                nge[j] <- .getResults(mdt, 1, "nge")
                pval[j] <- .getResults(mdt, 1, "pval")
                mcoef[j] <- paste(.getResults(mdt, 1, "ix"), collapse=",", sep="")
              } else {
                nge[j] <- getNExceeding(pmx, cvars[1])
                pval[j] <- getPValue(pmx, cvars[1])
                mcoef[j] <- cvars[1]
              }
              ##############################################
            }

            gmpmMainSum <- data.frame(mcoef, nge, pval, .getSig(pval))
            colnames(gmpmMainSum) <- c("Coef","N>=Orig","p-value", " ")
            if (nTests > 1) {
              rownames(gmpmMainSum) <- colnames(faclist)
            } else {
              rownames(gmpmMainSum) <- x@ivars
            }

            if (length(x@covars) > 0) {                                       
              vv <- faclist[!(rownames(faclist) %in% x@ivars),]
              if (!is.vector(vv)) {
                vv <- as.vector(colSums(vv))
                vv[vv>1] <- 1
              }
              gmpmMainSum <- gmpmMainSum[order(vv),]
            }

            #print("... exiting .mainSumProc ...")               
            return(gmpmMainSum)
          })

setMethod("getMainSummary",
          signature(x="GMPM"),
          function(x, byCovar=FALSE) {
            gg <- .prepareMainSum(x, byCovar)
            faclist <- gg[["faclist"]]
            allvars <- gg[["allvars"]]
            nTests <- gg[["nTests"]]
            #print("~~~ in getMainSummary (GMPM) ~~~")
            ff <- list(.mainSumProc(x, faclist, allvars, nTests, x@psec))
            names(ff) <- c("Main Results")
            #print("... exiting getMainSummary GMPM) ...")
            return(ff)
          })

setMethod("getMainSummary",
          signature(x="GMPM.mul"),
          function(x, byCovar=FALSE) {
            #print("~~~ in getMainSummary (GMPM.mul) ~~~")
            gg <- .prepareMainSum(x, byCovar)
            faclist <- gg[["faclist"]]
            allvars <- gg[["allvars"]]
            nTests <- gg[["nTests"]]
            mlist <- list()
            DVlevels <- .getDVlevels(x)
            nDVlevels <- dim(x@coef0)[1]            
            mnames <- paste(DVlevels[2:length(DVlevels)], DVlevels[1],
                  sep=" versus ")
            for (i in 1:nDVlevels) {
              psec <- .collapseMultinomPmx(x, i)
              mlist[[i]] <-
                .mainSumProc(x, faclist, allvars, nTests, psec)
            }
            names(mlist) <- mnames
            #print("... exiting getMainSummary (GMPM.mul) ...")
            return(mlist)
          })

setMethod(".getIXfromIV",
          signature(x="GMPM"),
          function(x, ivinc, includeCovars=TRUE) {
            flab <- x@IVcoef
            if (includeCovars) {
              if (length(x@covars) > 0) {
                for (i in 1:length(x@covars)) {
                  flab[[x@covars[i]]] <- x@covars[i]
                }
              } else {}
            } else {}
            
            cvars <- c()
            for (k in 1:length(x@coefTerms)) {
              ctmp <- c()
              for (m in 1:length(ivinc)) {
                ivForms <- flab[[ivinc[m]]]
                isIn <- ivForms %in% x@coefTerms[[k]]
                if (sum(isIn) > 0) {
                  if (sum(isIn) != 1) {
                    stop("something's wrong here; bailing out.")
                  } else {
                    ctmp <- c(ctmp, ivForms[isIn])
                  }
                }
              }
              if ((length(ctmp) == length(ivinc)) &&
                  (length(ctmp) == length(x@coefTerms[[k]]))) {
                cvars <- c(cvars, k)
              } else {}
            }
            return(cvars)
          })

           
setMethod(".getCoefTerms",
          signature(x="GMPM"),
          function(x, mycoef) {
            if (missing(mycoef)) {
              mycoef <- x@coef0
            }
            if (is.matrix(mycoef)) {
              n <- colnames(mycoef)
            } else {
              n <- names(mycoef)
            }
            nl <- list()
            for (i in 1:length(n)) {
              nl[[i]] <- strsplit(n[i], ":")[[1]]
            }
            x@coefTerms <- nl
            return(nl)
          }
          )

setMethod(".storeFitResult",
          signature(x="GMPM"),
          function(x, fit, section, index) {
            nameObject <- deparse(substitute(x))

            myFit <- coef(fit)
            x@psec[[section]][index,] <- myFit
                                        #x@pmx[index,] <- myFit

            assign(nameObject, x, envir=parent.frame())            
          }
          )

setMethod(".storeFitResult",
          signature(x="GMPM.mul"),
          function(x, fit, section, index) {
            
            nameObject <- deparse(substitute(x))

            myFit <- coef(fit)
            x@psec[[section]][index,,] <- myFit
            if (!is.null(fit$convergence)) {
              x@convergence <- if(fit$convergence==1) FALSE else TRUE
            }

            assign(nameObject, x, envir=parent.frame())            
          }
          )

setMethod(".reportProgress",
          signature(x="GMPM"),
          function(x, ix, maxruns, elapsed) {
            fmtstr1 <- paste("%",nchar(as.character(maxruns)),"d",
                            sep="")
            ixx <- ifelse( (ix*x@nCores)>x@gmpmControl[["maxruns"]],
                          x@gmpmControl[["maxruns"]], ix*x@nCores)
            fmtstr1.2 <- paste("%", nchar(as.character(x@gmpmControl[["maxruns"]])),"d", sep="")
            fmtstr2 <- paste(fmtstr1, "/", maxruns, " (", fmtstr1.2,
                             "/", x@gmpmControl[["maxruns"]],
                             ") ", ":", sep="")
            stmp1 <- sprintf(fmtstr2, ix, ixx, maxruns)
            #stmp2 <- sprintf("% 1.3f", tmp1)
            if (length(elapsed) > 1) {
              efac <- mean(elapsed)
            } else {
              efac <- elapsed[1]
            }
            cat(stmp1, sprintf("%1.3fs/sweep, ", efac))

            totsecs <- efac*(maxruns-ix)
            cat(sprintf("%02dh:", floor(totsecs / 3600)))
            if (totsecs >= 3600) {
              remn <- totsecs %% 3600
            } else {
              remn <- totsecs
            }
            cat(sprintf("%02dm:", floor(remn/60)))
            cat(sprintf("%02ds", round(remn %% 60)), "left\n")
            flush.console()
          }
          )


setMethod(".createPermMx",
          signature(x="GMPM"),
          function(x, maxruns) {
            # create matrix
            pmx <- 
              matrix(nrow=maxruns+1, ncol=length(x@coef0))
            pmx[1,] <- x@coef0
            colnames(pmx) <- names(x@coef0)
            x@pmx <- pmx

            x@psec <- .createMatrixSections(x, pmx)
            
            return(x@psec)
          }
          )

setMethod(".createPermMx",
          signature(x="GMPM.mul"),
          function(x, maxruns) {
            nameObject <- deparse(substitute(x))

            pmx <- array(dim=c(maxruns+1, dim(x@coef0)[1], dim(x@coef0)[2]))
            capture.output(f0 <- origFit(x))
            pmx[1,,] <-coef(f0)
              
            dimnames(pmx) <- list(run=1:(maxruns+1),
                                    dv=rownames(x@coef0),
                                    coef=colnames(x@coef0))
            x@pmx <- pmx

            if (!is.null(f0$convergence)) {
              x@convergence <- rep(FALSE, maxruns+1)
              x@convergence[1] <- if(f0$convergence == 1) FALSE else TRUE
            } else {
              warning("You are using an older version of package 'nnet' ",
                      "that doesn't provide\n", "information about ",
                      "convergence.  Please consider updating.\n",
                      "(see ?update.packages() for instructions).")
              x@convergence <- rep(TRUE, maxruns+1)              
            }

            x@psec <- .createMatrixSections(x, pmx)

            return(x@psec)
          }
          )

setMethod("getNExceeding",
          signature(x="GMPM"),
          function(x, y) {
            if (missing(y)) {
              ix <- .getIVix(x)
            } else {
              ix <- y
            }
            if (is.null(x@coef0) || is.null(x@pmx)) {
              stop("Model must be estimated (gmpmEstimate) before extracting results.")
            }
            if (length(setdiff(ix, x@ivix)) > 0) {
              print(names(x@coef0)[setdiff(ix, x@ivix)])
              stop("Error: p-values only meaningful for independent variables and IV*covariate interactions.")
            }

            return(getNExceeding(x@pmx, ix))
          }
          )

setMethod("getNExceeding",
          signature(x="matrix"),
          function(x, y) {
            pmx <- x
            coef0 <- x[1,]
            ctmp <- c()
            for (i in 1:length(y)) {
              ctmp <- c(ctmp,
                        sum(abs(pmx[,y[i]])>=abs(coef0[y[i]])))
            }
            
            return(ctmp)
          })


setMethod("getPValue",
          signature(x="GMPM"),
          function(x, y) {
            if (missing(y)) {
              y <- .getIVix(x)
            }
            return(getPValue(getPermMx(x), y))
            #nge <- getNExceeding(x, ix)
            #return(nge/(x@ncomp+1))
          }
          )

setMethod("getPValue",
          signature(x="matrix"),
          function(x, y) {
            nge <- getNExceeding(x, y)
            ntot <- dim(x)[1]
            return(nge/ntot)
          })

setMethod("gmpmCoef",
          signature(x="GMPM"),
          function(x) {
            if (is.null(x@coef0)) {
              x@coef0 <- coef(fitOnce(x))
            }
            return(x@coef0)
          }
          )

#setMethod("coef",
#          signature(x="GMPM"),
#          function(x) {
#            coefficients(x)
#          })

setMethod("getPermMx",
          signature(x="GMPM"),
          function(x) {return(x@psec)})

setMethod("permute",
    signature(x = "GMPM", thisiv = "character"),
    function (x, thisiv) 
          {
            if (x@nBetween > 1) {
              pa <- x@psBetween[["tiers"]][[thisiv]]
              xt1 <- x@psBetween[["freqs"]]
              nTiers <- length(pa)
              iv.p <- xt1[,thisiv]
              dfx <- x@df1

              nLevels <- length(levels(dfx[,thisiv]))
              rix <- sample(1:(nLevels*x@minN))
                                        #rix <- c(7:9,1:3,4:6)
              for (i in 1:nTiers) {
                thisTier <- pa[[i]]
                osel <- rep(NA,length(rix))
                for (j in 1:length(thisTier)) {
                  ix1 <- (j-1)*x@minN+1
                  ix2 <- ix1+x@minN-1
                  osel[ix1:ix2] <- sample(thisTier[[j]],x@minN,replace=F)
                }
                iv.p[osel] <- xt1[osel,thisiv][rix]
              }
              iv.p2 <- rep(iv.p, xt1[,ncol(xt1)])
              contrasts(iv.p2) <- contrasts(xt1[,thisiv])
              x@df1[,thisiv] <- iv.p2
            } else {
              if (x@nBetween == 1) {
                # only one between variable; this is easy
                frq <- x@psBetween[["freqs"]]
                iv.p0 <- sample(frq[,x@ivBetween[1]])
                iv.p <- rep(iv.p0,frq[,ncol(frq)])
                contrasts(iv.p) <- contrasts(x@df1[,x@ivBetween[1]])
                x@df1[,x@ivBetween[1]] <- iv.p
              } else {}
            }
            
            return(invisible(x))
          }
          )

setMethod("permute",
    signature(x = "GMPM", thisiv = "missing"),
    function (x) 
          {
            
            frq <- x@psWithin[["freqs"]]
            nReps <- x@psWithin[["nReps"]]
            for (i in 1:x@nWithin) {
              tiv <- x@ivWithin[i]
              tsch <- x@psWithin[["scheme"]][[tiv]]
              psch <- apply(tsch,2,sample)
              iv.p0 <- frq[,tiv]
              contrasts(iv.p0) <- contrasts(x@df1[,tiv])

              if (x@nBetween > 0) {
                psmx <- apply(x@psWithin[["smx"]],
                              2,permuteNA)[1:x@minN,]

                # for each between-s cell
                for (j in 1:ncol(psmx)) {
                  for (k in 1:ncol(psch)) {
                    sid <- x@psWithin[["cellfreqs"]][psmx[k,j],x@munit]
                    lvec <- frq[,x@munit]==sid
                    iv.p0[lvec] <- rep(psch[,k],each=nReps[tiv])
                  }
                }
              } else {
                iv.px <- as.vector(apply(psch,2,rep,x@nwrep[i]))
                iv.p0 <- rep(iv.px,each=nReps[tiv])
                iv.p0 <- factor(iv.p0, levels=levels(x@df1[,tiv]))
                contrasts(iv.p0) <- contrasts(x@df1[,tiv])
              }

              iv.p <- rep(iv.p0, frq[,ncol(frq)])
              contrasts(iv.p) <- contrasts(x@df1[,tiv])
              x@df1[,tiv] <- iv.p
            }

            return(invisible(x))
          }
          )

###############################

setMethod("origFit",
          signature(object = "GMPM"),
          function (object) {
            return(eval(as.call(object@fitcall)))
          }
          )

setMethod(".getIVix",
          signature(x = "GMPM"),
          function (x) {
            return(object@ivix)
          }
          )

          #####################

setMethod(".parseFormula",
          signature(object = "GMPM"),
          function(object, formula)
          {
            #print("~~~ in .parseFormula ~~~")
            nameObject <- deparse(substitute(object))

            fstr <- deparse(formula, width.cutoff=500L)
            if (length(fstr) > 1) {
              fstr <- paste(fstr, collapse="")
            } else {}
            object@mform <- formula

            fterms <- strsplit(fstr, "~")
            if (length(fterms[[1]]) != 2) {
              stop("Formula", paste("'", fstr, "'", sep=""),
                   "is malformed.  Must have exactly two sides, separated by '~'.")
            }
            lhs <- fterms[[1]][1]
            rhs <- fterms[[1]][2]
            spos <- gregexpr("\\|", rhs)[[1]][1]
            if (spos > 0) {
              rhs <- substr(rhs, 1, spos-1)
            }
            object@dform <- as.formula(paste(lhs, rhs, sep="~"))
            
            assign(nameObject, object, envir=parent.frame())
            return(invisible(object))
          }
          )

setMethod(".checkMultilevel",
          signature(object="GMPM"),
          function(object)
          {
            #print("~~~ in .checkMultilevel ~~~")            
            nameObject <- deparse(substitute(object))

            ## MULTILEVEL DATA ? ###################################
            object@nunits <- 1
            object@munit <- ""
            ff <- object@mform
            if ("|" %in% as.character(ff[[3]]))
              {
                object@munit <- as.character(ff[[3]])[3]
                object@nunits <- dim(xtabs(as.formula(
                                                      paste("~",
                                                            object@munit,
                                                            sep="")),
                                           data=object@df1))

                if (object@nunits == dim(object@df1)[1])
                  {
                    #cat("Warning: Only one observation per sampling unit; data are single-level.\nGrouping factor will be ignored.", "\n")
                    #stop("only one observation per sampling unit; your data are single-level.\ngmpm not yet configured for single-level data; please wait until next version.")
                    #object@nunits <- 1
                    #object@mform <- object@dform
                    #object@munit <- ""
                  } else {}
              } else {
                stop("No grouping factor supplied.\n gmpm not configured to handle single-level data; please wait for next version.")
                cat("No conditioning variable specified in model; assuming data are not multilevel.\n")
              }

            assign(nameObject, object, envir=parent.frame())
            return(invisible(object))
          }
          )

setMethod(".getDesign",
          signature(object="GMPM"),
          function(object)
          {
            #print("~~~ in .getDesign ~~~")
            nameObject <- deparse(substitute(object))
            
            mfactors <- attr(attr(object@df1, "terms"), "factors")
            IVnames <- object@ivars
            covars <- setdiff(rownames(mfactors)[-1], IVnames)
            if (length(covars) > 0) {
              object@covars <- covars
            }
            
            if (length(IVnames)==0)
              {
                cat("Warning: No independent variables specified;\n",
                    "will perform one-sample test.\n")
              }
            else
              {
                allfacts <- rownames(mfactors)[2:(dim(mfactors))[1]]
                if (length(intersect(IVnames, allfacts)) < length(IVnames)) {
                  stop("Malformed formula: one or more IVs in list",
                       "missing from formula.")
                }

                ## check whether IVs are all coded as factors
                IVinfo <- data.frame(nLevels=rep(NA, length(IVnames)),
                                     Type=rep(NA, length(IVnames)),
                                     Levels=rep(NA, length(IVnames)))
                rownames(IVinfo) <- IVnames


                for (i in 1:length(IVnames)) {
                  f1 <- object@df1[,IVnames[i]]
                  if (!is.factor(f1)) {
                    errstr <- paste("Warning: Converting ", IVnames[i],
                                    " to a factor")
                    cat(errstr, "\n")
                    object@df1[,IVnames[i]] <- factor(f1)
                  }
                  nLevels <- length(levels(object@df1[,IVnames[i]]))
                  IVinfo[IVnames[i],"nLevels"] <- nLevels
                  if (nLevels == 2)
                    {
                      cmx <- c(-.5,.5)
                    }
                  else
                    {
                      cmx <- matrix(nrow=nLevels, ncol=nLevels)
                      cmx[] <- 0
                      diag(cmx) <- 1
                      cmx <- cmx[,2:nLevels]
                      for (j in 2:nLevels) {
                        cmx[,(j-1)] <- cmx[,(j-1)]-mean(cmx[,(j-1)])
                      }
                    }
                  object@df1[,IVnames[i]] <- C(object@df1[,IVnames[i]], cmx)

                  IVinfo[IVnames[i],"Levels"] <-
                    paste(levels(object@df1[,IVnames[i]]),
                          collapse=",")

                }


                IVinfo$Type <- .getIVtypes(object@df1, object@ivars,
                                          object@munit)
                object@nWithin <- sum(IVinfo$Type=="within")
                object@nBetween <- sum(IVinfo$Type=="between")                
                object@IVinfo <- IVinfo
              }

            assign(nameObject, object, envir=parent.frame())
            return(invisible(object))
          }          
          )

setMethod(".buildFitCall",
          signature(object="GMPM"),
          function(object, arg.exclude=c(), ocall)
          {
            #print("~~~ in .buildFitCall (GMPM)~~~")            
            # build function call (common to all classes) ################
            # this function is only invoked by the '.buildFitCall'
            # functions for the child classes (GMPM.mul, GMPM.user).
            #
            # The current function handles situations common to all
            # classes, and passes the result back to the child class.
            #
            # Note that this function does not store the call in the object;
            # this is handled by the 'child' class.

            gmpm.args <- c("gmpmControl", "ivars")
            # note: update the above line when new arguments are
            # introduced to GMPM to avoid them being passed along
            # to fitting function.
            
            # get rid of gmpm-specific arguments
            ocall <- ocall[c(TRUE, !(names(ocall[2:length(ocall)]) %in%
                               c(arg.exclude, gmpm.args)))]
            mcl <- c(list(fun="user"), as.list(ocall[2:length(ocall)]))
            mcl$formula <- object@dform
            mcl$data <- quote(object@df1)
            #print(names(mcl))

            return(invisible(mcl))
          }
          )

setMethod(".buildFitCall",
          signature(object="GMPM.glm"),
          function(object, arg.exclude=c(), ocall)
          {
            #print("~~~ in .buildFitCall (GMPM.glm)~~~")                        

            ncall <- callNextMethod(object=object, arg.exclude=c(),
                           ocall=ocall)
            ncall[[1]] <- glm
            return(invisible(ncall))
          }            
          )

setMethod(".buildFitCall",
          signature(object="GMPM.mul"),
          function(object, arg.exclude=c(), ocall)
          {
            #print("~~~ in .buildFitCall (GMPM.mul)~~~")
            # build function call (common to all classes) ################

            ncall <- callNextMethod(object=object, arg.exclude=c("family"),
                           ocall=ocall)
            library(nnet)
            ncall[[1]] <- multinom
            return(invisible(ncall))
          }            
          )

setMethod(".buildFitCall",
          signature(object="GMPM.user"),
          function(object, arg.exclude=c(), ocall)
          {
            #print("~~~ in .buildFitCall (GMPM.user)~~~")

            ncall <- callNextMethod(object=object, arg.exclude=c("family"),
                           ocall=ocall)
            
            stop("user-defined GMPM objects not yet implemented.\n",
                 "Please wait for next version of gmpm package.")
          }            
          )

setMethod("fitOnce",
          signature(object="GMPM"),
          function(object)
          {
                                        ##print("~~~ in fitOnce (GMPM) ~~~")
            ff <- eval(as.call(object@fitcall))
            return(ff)
          }
          )

setMethod("fitOnce",
          signature(object="GMPM.mul"),
          function(object)
          {
                                        ##print("~~~ in fitOnce (GMPM) ~~~")
            capture.output(fit1 <- eval(as.call(object@fitcall)))
            if (!is.null(fit1$convergence)) {
              if (fit1$convergence == 1) {
                cat("Warning: multinom function did not converge.\n",
                    "Try increasing number of multinom iterations",
                    " (use 'maxit' argument).\n")
              } else {}
            } else {}
                
            return(fit1)
          }
          )

setMethod(".preparePermScheme",
          signature(object="GMPM"),
          function(object)
          {
            #print("~~~ in .preparePermScheme ~~~")
            nameObject <- deparse(substitute(object))            
            
            # prepare permutation scheme ####################
            # and re-sort the data set.

            # save typing/make code more readable
            nWithin <- object@nWithin
            nBetween <- object@nBetween
            rnames <- rownames(object@IVinfo)
            ivWithin <- rnames[object@IVinfo$Type=="within"]
            ivBetween <- rnames[object@IVinfo$Type=="between"]
            object@ivWithin <- ivWithin
            object@ivBetween <- ivBetween
            object@psBetween <- list()
            object@psWithin <- list()
            
            id <- object@munit
            x <- object@df1
            
            # sort data by within-subjects variables
            # (this makes perm tests cleaner)            
            x <- sortDF(x, c(ivBetween, id, ivWithin))
            object@df1 <- x

            # create permutation scheme for between subject variables
            # create summary table (xt1) for use in rearranging labels
            if (nBetween > 0) {
              xt1 <- combFreqs(x, c(ivBetween, id))
              xt2 <- combFreqs(xt1, ivBetween)
              object@minN <- minN <- min(xt2$N)
              object@psBetween[["freqs"]] <- xt1
              
              if (nBetween > 1) {
                pTiersB <- list()
                for (i in 1:nBetween) {
                  lvl <- levels(xt2[,ivBetween[i]])
                  xt3 <- combFreqs(xt1, setdiff(ivBetween, ivBetween[i]))
                  tlist <- list()
                  for (j in 1:nrow(xt3)) {

                                        # create tier name
                    flvl <- xt3[j,1:(ncol(xt3)-1)]
                    labs <- flvl
                                        # defactorize

                    if (is.data.frame(labs)) {
                      for (p in 1:ncol(labs)) { # error when labs is vector
                        labs[,p] <- levels(labs[,p])[labs[,p]]
                      }
                    } else {
                      if (is.factor(labs)) {
                        labs <- levels(labs)[labs]
                        names(labs) <- colnames(xt3)[1]
                      } else {}
                    }
                    tname <- paste(labs,collapse=":")                    

                    llist <- list()
                    for (k in 1:length(lvl)) {
                      flist <- list(lvl[k])
                      names(flist) <- ivBetween[i]
                      flist <- c(flist,as.list(labs))
                      llist[[lvl[k]]] <- getDFix(xt1,flist)
                    }
                    tlist[[tname]] <- llist
                  }
                  pTiersB[[ivBetween[i]]] <- tlist
                }
                object@psBetween[["tiers"]] <- pTiersB
              } else {}
              nSubjRand <- minN
            } else {
              nSubjRand <- object@nunits
            }

            # generate permutation scheme for within factors
            # TO DO: check the design for missing factor levels
            # within subject (could be a problem)
            if (object@nWithin > 0) {
              xt4 <- combFreqs(x, c(ivBetween, id))
              xt5 <- combFreqs(x, c(ivBetween, id, ivWithin))
              wlist <- list()
              nReps <- rep(1,nWithin)
              names(nReps) <- ivWithin
              nwrep <- rep(1, nWithin)
              nwLevels <- rep(NA, nWithin)
              for (i in 1:nWithin) {
                tiv <- ivWithin[i]
                nwLevels[i] <- nLevels <- length(levels(x[,tiv]))
                
                if (i == 1) {
                  nwrep[i] <- 1
                } else {
                  for (j in 1:(i-1)) {
                    nwrep[i] <- nwrep[i]*nwLevels[j]
                  }
                }
                
                mx <- matrix(levels(x[,tiv]),nrow=nLevels,ncol=nSubjRand)
                wlist[[ivWithin[i]]] <- mx
                if (i < nWithin) {
                  for (j in (i+1):nWithin) {
                    nReps[i] <- nReps[i]*length(levels(x[,ivWithin[j]]))
                  }
                } else {}
              }
              object@nwrep <- nwrep
              object@psWithin$scheme <- wlist
              object@psWithin$freqs <- xt5
              object@psWithin$nReps <- nReps
              
              if (nBetween > 0) {
                maxN <- max(xt2$N)
                smx <- matrix(nrow=maxN,ncol=nrow(xt2))
                for (i in 1:nrow(xt2)) {
                  ix <- xt2[i,-ncol(xt2)]
                  if (is.factor(ix)) {
                    ix <- levels(ix)[ix]
                    names(ix) <- colnames(xt2)[-ncol(xt2)]
                  } else {}
                  dfix <- getDFix(xt1,ix)
                  smx[1:length(dfix),i] <- dfix
                }
                object@psWithin[["cellfreqs"]] <- xt1
                object@psWithin[["smx"]] <- smx
              } else {}
            }

            assign(nameObject, object, envir=parent.frame())
            return(invisible(object))            
          }
          )

gmpFit <- function(x,y=NULL) {
  stop("Function gmpFit superseded by gmpmEstimate as of version 0.4-0.")
}

setMethod(".getIVix",
          signature(x="GMPM"),
          function(x) {
            ivix <- c()
            n <- names(x@coef0)
            nl <- x@coefTerms

            ivix <- c()
            for (i in 1:length(nl)) {
              varsIn <- nl[[i]]
              for (j in 1:length(x@IVcoef)) {
                for (m in 1:length(x@IVcoef[[j]])) {
                  if (x@IVcoef[[j]][m] %in% varsIn) {
                    ivix <- c(ivix, i)
                  }
                }
              }
            }
            ivix <- sort(unique(ivix))
            
            return(ivix)
          }
          )


setMethod("permSpace",
          signature(object="GMPM"),
          function(object)  {
            psBetween <- object@psBetween
            nCellsPerUnit <- object@nCellsPerUnit
            ivBetween <- object@ivBetween
            ivWithin <- object@ivWithin
            nunits <- object@nunits
            recChoose <- function(v) {
              tot <- sum(v)
              if (length(v)==1) {
                return(1)
              } else {
                n1 <- choose(tot, v[1])
                return(n1*recChoose(v[2:length(v)]))
              }
            }
            
            if (length(ivBetween)) {
              betform <- as.formula(paste("Freq ~", paste(ivBetween,collapse="+",sep="")))
              xt1 <- as.data.frame(xtabs(betform, data=psBetween))
              xt1 <- xt1[xt1$Freq > 0,] # get rid of unused levels
              bfac <- recChoose(xt1$Freq)
            } else {
              bfac <- 1
            }

            if (length(ivWithin)) {
              cc <- 1
              for (j in 1:length(ivWithin)) {
                cc <- cc*length(levels(object@df1[,ivWithin[j]]))
              }
              wfac <- cc^nunits
            } else {
              wfac <- 1
            }

            return(bfac*wfac)
          }
          )

setMethod("coefNames",
          signature(x="GMPM"),
          function(x)  {
            return(names(x@coef0))
          })

setMethod("coefNames",
          signature(x="GMPM.mul"),
          function(x)  {
            return(colnames(x@coef0))
          })

setMethod("gmpmEstimate",
          signature(x="GMPM"),
          function(x,gmpmControl)
          {
            #print("~~~ in gmpmFit (GMPM) ~~~")
            if (!missing(gmpmControl)) {
              .setOpts(x, gmpmControl)
              #gmpmControl <- x@gmpmControl
            } else {}

            x@nCores <- .calculateCores(x@gmpmControl)

            if (is.null(x@coef0)) {
              x@coef0 <- coef(fitOnce(x))
            }

            maxruns <- x@gmpmControl[["maxruns"]]
            report.interval <- x@gmpmControl[["report.interval"]]
            outfile <- x@gmpmControl[["outfile"]]
            
            elapsed <- rep(NA, maxruns)

            if (!is.null(outfile)) {
              stop("writing to outfile not yet supported; bailing out.\n")
              .writeFit(x, coef(fitOnce(x)), outfile, FALSE)
            } else {}

            # calculate number of estimation sections
            x@psec <- .createPermMx(x, maxruns)            
            x@nSections <- x@nBetween + (x@nWithin > 0)
            if (x@nBetween > 0) {
              pbnames <- names(x@psec)[1:x@nBetween]
            } else {}

            # one "sweep" is defined as a set of parallel runs
            # on each of the requested processing cores.
            #
            # the following code sets up the calls for doing the
            # fitting over multiple processors.
            nsweeps <- floor(maxruns / x@nCores)
            listsize <- rep(x@nCores, nsweeps)
            if ((maxruns %% x@nCores) > 0) {
              nsweeps <- nsweeps + 1
              listsize <- c(listsize, maxruns %% x@nCores)
            } else {}

            if ("parallel" %in% installed.packages()) {
              mycall <- "mclapply"
            } else {
              mycall <- "lapply"
            }
            
            x.list <- list()
            for (i in 1:x@nCores) {
              x.list[[i]] <- x
            }
            
            for (i in 1:nsweeps) {
              if (i == nsweeps) {
                x.list <- x.list[1:listsize[i]]
              } else {}
              thisrow <- (i-1)*x@nCores+2
              
              t1 <- proc.time()["elapsed"]
              
              if (x@nBetween > 0) {
                for (j in 1:x@nBetween) {
                  x2.list <- do.call(mycall, args=list(X=x.list, FUN=permute, thisiv=x@ivBetween[j]))
                  myFit.list <- do.call(mycall, args=list(X=x2.list, FUN=fitOnce))
                  for (mm in 1:listsize[i]) {
                    .storeFitResult(x, myFit.list[[mm]], pbnames[j], thisrow+mm-1)
                  }
                  #myFit <- fitOnce(x2)
                  #x@psec[[pbnames[j]]][i+1,] <- cfs
                }
              } else {}

              if (x@nWithin > 0) {
                x2.list <- do.call(mycall, args=list(X=x.list, FUN=permute))
                #x2 <- permute(x)
                myFit.list <- do.call(mycall, args=list(X=x2.list, FUN=fitOnce))
                for (mm in 1:listsize[i]) {
                  .storeFitResult(x, myFit.list[[mm]], "within", thisrow+mm-1)
                }
                #x@psec[["within"]][i+1,] <- cfs
                #.storeFitResult(x, ff, i+1)
              } else {}

              #if (!is.null(outfile)) {
              #  .writeFit(x, myFit, outfile, TRUE)
              #} else {}

              t2 <- proc.time()["elapsed"]
              elapsed[i] <- t2-t1

              if (report.interval > 0) {
                if ((i == maxruns) ||
                    ((i %% report.interval)==0)) {
                  .reportProgress(x, i, nsweeps, elapsed[1:i])
                } else {}
              } else {}
            }

            x@ncomp <- sum(listsize[1:i])
            x@ndigits <- ifelse(x@ncomp > 9, ceiling(log(x@ncomp+1,base=10)), 1)

            #print("~~~ leaving gmpmFit ~~~")            
            return(x)
                        
          })

setMethod(".createMatrixSections",
          signature(x="GMPM"),
          function(x,pmx) {
            psec <- list()
            # now create matrix sections
            x@nSections <- x@nBetween + (x@nWithin > 0)
            sb <- list()
            if (x@nBetween > 0) {
              for (i in 1:x@nBetween) {
                sb[[i]] <- pmx
              }
              names(sb) <-
                paste("between",rownames(x@IVinfo)[x@IVinfo$Type=="between"],
                      sep=":")
              psec <- c(sb)
            }
                
            if (x@nWithin > 0) {
              psec[["within"]] <- pmx
            } else {}

            x@psec <- psec

            return(psec)            
          })

setMethod(".collapseMultinomPmx",
          signature(x="GMPM.mul"),
          function(x,index) {
              psec <- list()
              for (j in 1:length(x@psec)) {
                psec[[j]] <- x@psec[[j]][,index,]
              }
              names(psec) <- names(x@psec)
              return(psec)
          })

setMethod(".setOpts",
          signature(x="GMPM"),
          function(x, opts) {
            nameObject <- deparse(substitute(x))            
            
            if (!is.null(opts[["nCores"]])) {
              x@gmpmControl[["nCores"]] <- .calculateCores(opts)
              x@nCores <- as.integer(x@gmpmControl[["nCores"]])
            } else {}

            if (!is.null(opts[["maxruns"]])) {
              x@gmpmControl[["maxruns"]] <- opts[["maxruns"]]
            } else {}

            if (!is.null(opts[["report.interval"]])) {
              x@gmpmControl[["report.interval"]] <- opts[["report.interval"]]
            } else {}

            if (!is.null(opts[["outfile"]])) {
              x@gmpmControl[["outfile"]] <- opts[["outfile"]]
            } else {}
            assign(nameObject, x, envir=parent.frame())            
            return(invisible(x))
          })
