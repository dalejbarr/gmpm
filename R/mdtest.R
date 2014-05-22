setClass("Mdtest",
         representation(
                        nSections="numeric", # number of sections
                        mdsections="list" # contains all results
                        ),
         prototype(nSections=1)
         )

setClass("Mdtest.sum",
         representation(
                        )
         )

setGeneric(
           name="mdTest",
           def=function(x, y){standardGeneric("mdTest")}
           )

setGeneric(
           name=".getResults",
           def=function(x, ix, result, nSection=1){
             standardGeneric(".getResults")}
           )

setMethod("show",
          signature(object = "Mdtest"),
          function (object) 
          {
            summary(object)
          }
          )

setMethod("summary",
          signature(object = "Mdtest"),
          function (object, ...) {

            cat("Mahalanobis distance test", "\n\n")

            secnames <- names(object@mdsections)
            for (m in 1:object@nSections) {
              if (m > 1) {
                cat("\n")
              }
              cat(secnames[m], ":\n")
              mdresults <- object@mdsections[[m]]

              mx <- matrix(nrow=length(mdresults), ncol=3)
              rnames <- 1:length(mdresults)
              cnames <- rep(NA, length(mdresults))
              for (i in 1:length(mdresults)) {
                x <- mdresults[[i]]
                mx[i,] <- round(c(x$mdOrig,
                                  x$nge,
                                  x$pval),4)
                if (!is.null(attr(mdresults[[i]], "test"))) {
                  rnames[i] <- attr(mdresults[[i]], "test")
                } else {
                  rnames[i] <- ""
                }
                                        # enter coefficient names.
                #cnames[i] <- x[["cnames"]]
                cnames[i] <- paste(x$ix, collapse=",", sep="")
              }
              df1 <- data.frame(cnames, mx[,1], mx[,2], mx[,3], .getSig(mx[,3]))
              colnames(df1) <- c("Coef", "Distance", "N>=Orig", "p-value"," ")
              rownames(df1) <- rnames
              print(df1)
            }
          }
          )

setMethod("mdTest",
          signature(x = "GMPM", y="vector"),
          function(x, y) {
            stop("Sorry... this function not yet implemented.  Use the matrix version instead.")
            index <- y
            ###print("~~~ in mdTest (GMPM, vector) ~~~")            
            xmdt <- new("Mdtest")
            xmdt@nSections <- 1
            xmdt@mdsections[["Main Results"]] <-
              .mdTestCalc(getPermMx(x), index)
            return(xmdt)
          }
          )

setMethod("mdTest",
          signature(x = "GMPM.mul", y="vector"),
          function(x, y) {
            index <- y
                                        ##print("~~~ in mdTest (GMPM.mul, vector) ~~~")
            xmdt <- new("Mdtest")
            DVlevels <- levels(x@df1[,x@DVname])
            nDVlevels <- length(DVlevels)
            xmdt@nSections <- nDVlevels-1
            for (i in 1:xmdt@nSections) {
              sTitle <- paste(dimnames(x@pmx)[["dv"]][i],
                              "versus", DVlevels[1], sep=" ")
              xmdt@mdsections[[sTitle]] <-
                .mdTestCalc(x@pmx[,i,], index)
            }
                                        ##print("... exiting mdTest (GMPM.mul, vector) ~~~")
            return(xmdt)
          }
          )

setMethod("mdTest",
          signature(x = "matrix"),
          function(x, y) {
            index <- y
            xmdt <- new("Mdtest")
            xmdt@nSections <- 1
            xmdt@mdsections[["Main Results"]] <-
              .mdTestCalc(x, index)
            return(xmdt)
          })

setMethod(".getResults",
          signature(x = "Mdtest"),
          function(x, ix, result, nSection=1) {
            return(x@mdsections[[nSection]][[ix]][[result]])
          }
          )
                 

.mdTestCalc <- function(pmx0, index) {

            # internal function to perform mdTests
            # first parse the arguments
            # will return a list with results
            
            if (missing(index)) {
              stop("Argument 'index' must be supplied.\nSee (?mdTest) for details.")
            } else {
              if (is.character(index)) {
                # search first among IVs
                sstr <- index
                index <- match(sstr, colnames(pmx0), 0)
                index <- index[index > 0]
                mlist <- list(index)
                # later: search among factors
              } else {
                if (is.list(index)) {
                  mlist <- index
                } else {
                  if (!is.numeric(index)) {
                    stop("Argument 'index' not recognized.\n")
                  } else {
                    mlist <- list(index)
                  }
                }
              }
            }

            mdr <- list()
            for (j in 1:length(mlist)) {
              index <- mlist[[j]]
              if (is.character(index)) {
                sstr <- index
                cfnames <- colnames(pmx0)
                if (sum(sstr %in% cfnames) < length(sstr)) {
                  print(cfnames)
                  stop("Names '",
                       paste(sstr[!(sstr %in% cfnames)], collapse=", "),
                       "' from index list element #", j,
                       " were not found in original coefficients.")
                }
                index <- match(sstr, colnames(pmx0))
              } else {
              }
              
              mdr[[j]] <- list()
              mdr[[j]]$gmpmcoef <- pmx0[1,index]
              mdr[[j]]$ix <- index

              mdr[[j]]["cnames"] <- paste("(",
                                 paste(names(mdr[[j]]$gmpmcoef),
                                       collapse=", ",sep=""),
                                 ")", sep="")
              
              if (!is.null(names(mlist))) {
                attr(mdr[[j]], "test") <-
                  names(mlist)[j]
              }
              
              pmx <- pmx0[-1,index]
              if (length(index)==1) {
                mdr[[j]]$nruns <- length(pmx)
                x.mus <- mean(pmx)
                mdr[[j]]$mdOrig <- abs(mdr[[j]]$gmpmcoef - x.mus)
                mdr[[j]]$nge <-
                  sum(abs(pmx) >= abs(mdr[[j]]$gmpmcoef))
                mdr[[j]]$pval <-
                  mdr[[j]]$nge / (mdr[[j]]$nruns+1)
              } else {
                x.mus <- colMeans(pmx)
                x.cov <- cov(pmx)
                x.cov.inv <- solve(x.cov)
                mdr[[j]]$vcov <- x.cov.inv

                mdr[[j]]$mdOrig <-
                  as.numeric(sqrt(t(mdr[[j]]$gmpmcoef - x.mus) %*%
                                  x.cov.inv %*%
                                  (mdr[[j]]$gmpmcoef - x.mus)))
                mdr[[j]]$nruns <- length(pmx[,1])

                mdr[[j]]$mdPerm <- rep(0, mdr[[j]]$nruns+1)
                mdr[[j]]$mdPerm[1] <- mdr[[j]]$mdOrig
                for (i in 2:(mdr[[j]]$nruns+1)) {
                  mdr[[j]]$mdPerm[i] <-
                    sqrt(t(pmx[i-1,] - x.mus) %*%
                         x.cov.inv %*%
                         (pmx[i-1,] - x.mus))
                }
                mdr[[j]]$nge <-
                  sum(mdr[[j]]$mdPerm >= mdr[[j]]$mdOrig)
                mdr[[j]]$pval <-
                  mdr[[j]]$nge / (mdr[[j]]$nruns+1)
              }
            }

            return(mdr)
          }

