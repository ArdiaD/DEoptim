DEoptim.control <- function(VTR = -Inf, strategy = 2, bs = FALSE, NP =
                            NA, itermax = 200, CR = 0.5, F = 0.8,
                            trace = TRUE, initialpop = NULL,
                            storepopfrom = itermax + 1, storepopfreq =
                            1, p = 0.2, c = 0, reltol, steptol,
                            parallelType = c("none", "auto", "parallel", "foreach"),
                            cluster = NULL, packages = c(),
			    parVar = c(), foreachArgs = list(),
			    parallelArgs = NULL) {
  if (itermax <= 0) {
    warning("'itermax' <= 0; set to default value 200\n", immediate. = TRUE)
    itermax <- 200
  }
  if (F < 0 | F > 2) {
    warning("'F' not in [0,2]; set to default value 0.8\n", immediate. = TRUE)
    F <- 0.8
  }
  if (CR < 0 | CR > 1) {
    warning("'CR' not in [0,1]; set to default value 0.5\n", immediate. = TRUE)
    CR <- 0.5
  }
  if (strategy < 1 | strategy > 6) {
    warning("'strategy' not in {1,...,6}; set to default value 2\n",
            immediate. = TRUE)
    strategy <- 2
  }
    
  bs <- (bs > 0)

  if ( trace < 0 ) {
    warning("'trace' cannot be negative; set to 'TRUE'")
    trace <- TRUE
  }

  storepopfreq <- floor(storepopfreq)
  if (storepopfreq > itermax)
    storepopfreq <- 1

  if (p <= 0 || p > 1) {
    warning("'p' not in (0,1]; set to default value 0.2\n", immediate. = TRUE)
    p <- 0.2
  }

  if (c < 0 || c > 1) {
    warning("'c' not in [0,1]; set to default value 0\n", immediate. = TRUE)
    c <- 0
  }

  if (missing(reltol)) {
    reltol <- sqrt(.Machine$double.eps)
  }

  if (missing(steptol)) {
    steptol <- itermax
  }

  if(!(is.null(initialpop))) {
    if(is.na(NP))
      if(is.matrix(initialpop))
        NP <- dim(initialpop)[1]
      else stop("initialpop must be a matrix")
    else
      if(NP != dim(initialpop)[1]) {
        warning("Resetting NP to the number of rows in initialpop")
        NP <- dim(initialpop)[1]
      }
  }
 #####################
  # handle parallel options

  #check for a single parallelType
  if(missing(parallelType) || length(parallelType)>1){
      parallelType<-parallelType[1]
  }
  # handle 'auto' auto-detect
  if(parallelType=='auto'){
     pkgs<-.packages()
     rv<-R.Version()
     if('foreach' %in% pkgs){
         parallelType='foreach'
     } else if (('parallel' %in% pkgs) ||
                (rv$major>=2 && rv$minor>=14.2)
                ){
         parallelType='parallel'
     } else {
         parallelType='none'
     }
  }
  #support old deprecated parallelType arguments
  if(is.numeric(parallelType)) {
    parallelType <- switch(parallelType+1, 'none', 'parallel', 'foreach')
  }

  #handle deptrecated parallel arguments, set sensible defaults, etc.
  switch(parallelType,
          foreach = {
              if(missing(parallelArgs) && hasArg(foreachArgs)){
                  parallelArgs<-match.call(expand.dots=TRUE)$foreachArgs
              }
              if(is.null(parallelArgs$.packages) ){
                  if(hasArg(packages)) parallelArgs$.packages<-match.call(expand.dots=TRUE)$packages
              }
          },
          parallel = {
              if(missing(packages) || !hasArg(packages)){
                  packages<-(.packages())
              }
          }
  )
  # end parallel options
  ######################

  # format and return
  list(VTR = VTR, strategy = strategy, NP = NP, itermax = itermax, CR
       = CR, F = F, bs = bs, trace = trace, initialpop = initialpop,
       storepopfrom = storepopfrom, storepopfreq = storepopfreq, p =
       p, c = c, reltol = reltol, steptol = steptol, parallelType =
           parallelType, cluster = cluster,
       packages = packages, parVar = parVar, foreachArgs =
           foreachArgs, parallelArgs = parallelArgs)
}

DEoptim <- function(fn, lower, upper, control = DEoptim.control(), ...,
                    fnMap=NULL) {
    if (length(lower) != length(upper))
        stop("'lower' and 'upper' are not of same length")
    if (!is.vector(lower))
        lower <- as.vector(lower)
    if (!is.vector(upper))
        upper <- as.vector(upper)
    if (any(lower > upper))
        stop("'lower' > 'upper'")
    if (any(lower == "Inf"))
        warning("you set a component of 'lower' to 'Inf'. May imply 'NaN' results", immediate. = TRUE)
    if (any(lower == "-Inf"))
        warning("you set a component of 'lower' to '-Inf'. May imply 'NaN' results", immediate. = TRUE)
    if (any(upper == "Inf"))
        warning("you set a component of 'upper' to 'Inf'. May imply 'NaN' results", immediate. = TRUE)
    if (any(upper == "-Inf"))
        warning("you set a component of 'upper' to '-Inf'. May imply 'NaN' results", immediate. = TRUE)
    if (!is.null(names(lower)))
        nam <- names(lower)
    else if (!is.null(names(upper)) & is.null(names(lower)))
        nam <- names(upper)
    else
        nam <- paste("par", 1:length(lower), sep = "")
    
    ctrl <- do.call(DEoptim.control, as.list(control))
    ctrl$npar <- length(lower)
    if(is.na(ctrl$NP))
        ctrl$NP <- 10*length(lower)
    if (ctrl$NP < 4) {
        warning("'NP' < 4; set to default value 10*length(lower)\n",
                immediate. = TRUE)
        ctrl$NP <- 10*length(lower)
    }
    if (ctrl$NP < 10*length(lower)) 
        warning("For many problems it is best to set 'NP' (in 'control') to be at least ten times the length of the parameter vector. \n", immediate. = TRUE)
    if (!is.null(ctrl$initialpop)) {
        ctrl$specinitialpop <- TRUE
        if(!identical(as.numeric(dim(ctrl$initialpop)), as.numeric(c(ctrl$NP, ctrl$npar))))
            stop("Initial population is not a matrix with dim. NP x length(upper).")
    }
    else {
        ctrl$specinitialpop <- FALSE
        ctrl$initialpop <- 0.0
    }
    ##
    ctrl$trace <- as.numeric(ctrl$trace)
    ctrl$specinitialpop <- as.numeric(ctrl$specinitialpop)
    ctrl$initialpop <- as.numeric(ctrl$initialpop)
    if(!is.null(ctrl$cluster)) { ## use provided cluster
        if(!inherits(ctrl$cluster, "cluster"))
            stop("cluster is not a 'cluster' class object")
        parallel::clusterExport(ctrl$cluster, ctrl$parVar)
        if(!is.null(ctrl$parallelArgs))
          fnPop <- function(`*params`, ...) {
            parallel::parApply(cl=ctrl$cluster,X=`*params`,MARGIN=1,FUN=fn,ctrl$parallelArgs,...)
        }
	      else
	        fnPop <- function(`*params`, ...) {
            parallel::parApply(cl=ctrl$cluster,X=`*params`,MARGIN=1,FUN=fn,...)
        }
    }
    else if(ctrl$parallelType == 'foreach') { ## use foreach
      use.foreach <- 'foreach' %in% installed.packages() 
      if(!use.foreach)
         stop("foreach package not available but parallelType set to 'foreach'")
        if(!foreach::getDoParRegistered()) {
            foreach::registerDoSEQ()
        }
        args <-  ctrl$foreachArgs
        fnPop <- function(`*params`, ...) {
            my_chunksize <- ceiling(NROW(`*params`)/foreach::getDoParWorkers())
            my_iter <- iterators::iter(`*params`,by="row",chunksize=my_chunksize)
            args$i <- my_iter
            if(is.null(args$.combine)) args$.combine <- c
            if (!is.null(args$.export))
                args$.export = c(args$.export, "fn")
            else
                args$.export = "fn"
            if (is.null(args$.errorhandling))
                args$.errorhandling = c('stop', 'remove', 'pass')
            if (is.null(args$.verbose))
                args$.verbose = FALSE
            if (is.null(args$.inorder))
                args$.inorder = TRUE
            if (is.null(args$.multicombine))
                args$.multicombine = FALSE
            foreach::"%dopar%"(do.call(foreach::foreach, args), apply(X=i,MARGIN=1,FUN=fn,...))
      
        }
    }
    else if(ctrl$parallelType == 'parallel'){ ## use parallel 
        cl <- parallel::makeCluster(parallelly::availableCores())
        packFn <- function(packages) {
            for(i in packages)
                library(i, character.only = TRUE)
        }
        parallel::clusterCall(cl, packFn, ctrl$packages)
	if(is.null(ctrl$parVar)) ctrl$parVar <- ls()
        parallel::clusterExport(cl=cl, varlist=ctrl$parVar, envir = environment())
        fnPop <- function(`*params`, ...) {
            parallel::parApply(cl=cl,X=`*params`,MARGIN=1,FUN=fn,...)
        }
    }
    else {  ## use regular for loop / apply
        fnPop <- function(`*params`, ...) {
            apply(X = `*params`, MARGIN = 1, FUN = fn, ...)
        }
    }
    
  ## Mapping function
  fnMapC <- NULL
  if(!is.null(fnMap)) {   
    fnMapC <- function(`*params`,...) {
            ## run mapping function
            mappedPop <- t(apply(X = `*params`, MARGIN = 1, FUN = fnMap, ...))
            if(all(dim(mappedPop) != dim(`*params`)))  ## check results
                stop("mapping function did not return an object with ",
                     "dim NP x length(upper).")
            dups <- duplicated(mappedPop)  ## check for duplicates
            np <- NCOL(mappedPop)
            tries <- 0
      while(tries < 5 && any(dups)) {
        ##print('dups!'); flush.console()
        nd <- sum(dups)
        ## generate new random population member
        newPop <- matrix(runif(nd*np),ncol=np)
        newPop <- rep(lower,each=nd) + newPop * rep(upper-lower,each=nd)
        ## replace duplicate with _mapped_ random member
        mappedPop[dups,] <- t(apply(newPop, MARGIN = 1, FUN = fnMap))
        dups <- duplicated(mappedPop)  ## re-check for duplicates
        tries <- tries + 1
      }
      if(tries==5)
        warning("Could not remove ",sum(dups)," duplicates from the mapped ",
                "population in 5 tries. Evaluating population with duplicates.",
                call.=FALSE, immediate.=TRUE)
      ## memcpy fails if mappedPop isn't double (need TYPEOF switch in C?)
      storage.mode(mappedPop) <- "double"
      mappedPop
    }
  }

  outC <- .Call("DEoptimC", lower, upper, fnPop, ctrl, new.env(), fnMapC, PACKAGE="DEoptim")

  if(ctrl$parallelType == "parallel")
    parallel::stopCluster(cl) 
  
  if (length(outC$storepop) > 0) {
    nstorepop <- floor((outC$iter - ctrl$storepopfrom) / ctrl$storepopfreq)
    storepop <- list()
    cnt <- 1
    for(i in 1:nstorepop) {
      idx <- cnt:((cnt - 1) + (ctrl$NP * ctrl$npar))
      storepop[[i]] <- matrix(outC$storepop[idx], nrow = ctrl$NP, ncol = ctrl$npar,
                         byrow = TRUE)
      cnt <- cnt + (ctrl$NP * ctrl$npar)
      dimnames(storepop[[i]]) <- list(1:ctrl$NP, nam)
    }
  }
  else {
    storepop = NULL
  }

  names(outC$bestmem) <- nam
  iter <- max(1,as.numeric(outC$iter))

  names(lower) <- names(upper) <- nam
  bestmemit <- matrix(outC$bestmemit[1:(iter * ctrl$npar)],
                      nrow = iter, ncol = ctrl$npar, byrow = TRUE)

  dimnames(bestmemit) <- list(1:iter, nam)
  storepop <- as.list(storepop)

  outR <- list(optim = list(
                 bestmem = outC$bestmem,
                 bestval = outC$bestval,
                 nfeval = outC$nfeval,
                 iter = outC$iter),
               member = list(
                 lower = lower,
                 upper = upper,
                 bestmemit = bestmemit,
                 bestvalit = outC$bestvalit,
                 pop = t(outC$pop),
                 storepop = storepop))
  
  attr(outR, "class") <- "DEoptim"
  return(outR)
}