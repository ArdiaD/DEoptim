"DEoptim" <-
function(FUN, lower, upper, control = list(), ...) {
  ## Differential Evolution Optimization
  ## David Ardia -- 2005-08-14

  if (missing(FUN))
    stop("'FUN' is missing") 
  FUN <- match.fun(FUN)
  
  if (missing(lower) || missing(upper))
    stop("'lower' or 'upper' is missing")
  if (length(lower) != length(upper))
    stop("'lower' and 'upper' are not of same length")
  if (!is.vector(lower))
    lower <- as.vector(lower)
  if (!is.vector(upper))
    XVmax <- as.vector(upper)
  if (any(lower > upper))
    stop("'lower' > 'upper'")
  
  d <- length(lower)
  con <- list(VTR = 1.e-6, itermax = 200,
              NP = 50, F = 0.8, CR = 0.5, strategy = 7,
              refresh = 10)
  con[names(control)] <- control
  
  if (con$itermax <= 0) {
    warning("'itermax' <= 0; set to default value 200\n")
    con$itermax <- 200
  }

  if (con$NP < 1) {
    warning("'NP' < 1; set to default value 50\n")
    con$NP <- 50
  }
  NP = con$NP
  
  if (con$F < 0 | con$F > 2) {
    warning("'F' not in [0,2]; set to default value 0.8\n")
    con$F <- 0.8
  }

  if (con$CR < 0 | con$CR > 1) {
    warning("'CR' not in [0,1]; set to default value 0.5\n")
    con$CR <- 0.5
  }
  
  if (con$strategy < 1 | con$strategy > 10) {
    warning("'strategy' not in {1,...,10}; set to default value 7\n")
    con$strategy <- 7
  }
  
  con$refresh <- floor(con$refresh)
  if (con$refresh > con$itermax)
    con$refresh <- 0
  
  ## Sub-functions
  fn.zeros <- function(nr, nc)
    matrix(rep.int(0, nr * nc), nrow = nr)
  
  fn.checkBoundaries <- function(x, lower, upper) {
    r <- apply(rbind(lower, x), 2, max)
    r <- apply(rbind(upper, r), 2, min)
    r 
  }
  
  ##
  ## Initialize population and some arrays
  ##

  pop = matrix(rep.int(lower, NP), nrow = NP, byrow = TRUE) +
    matrix(runif(NP * d), nrow = NP) *
      matrix(rep.int(upper - lower, NP), nrow = NP, byrow = TRUE)
  
  popold <- fn.zeros(NP,d) ## toggle population
  val <- rep.int(0,NP) ## create and reset the "cost array"
  bestmem <- bestmemit <- rep.int(0,d) ## best population member ever and iteration

  ##
  ## Evaluate the best member after initialization
  ##

  nfeval <- NP ## number of function evaluations
  val <- apply(pop, 1, FUN, ...)
  bestval <- bestvalit <- min(val)
  ibest <- match(bestvalit, val)
  bestmem <- bestmemit <- pop[ibest,]
  
  ##
  ## DE - optimization
  ##
  ## popold is the population which has to compete. It is
  ## static through one iteration. pop is the newly emerging population.
  ##
  
  pm1 <- pm2 <- pm3 <- pm4 <- pm5 <- fn.zeros(NP,d) ## initialize population matrix 1 - 5
  bm <- ui <- mui <- mpo <- fn.zeros(NP,d)
  rot <- seq(from = 0, by = 1, to = (NP-1))## rotating index array (size NP)
  rotd <- seq(from = 0, by = 1, to = (d-1)) ## rotating index array (size d)
  rt <- fn.zeros(NP,NP) ## another rotating index array
  rtd <- fn.zeros(d,d) ## rotating index array for exponential crossover
  a1 <- a2 <- a3 <- a4 <- a5 <- fn.zeros(NP,NP) ## index array 1 - 5
  ind <- fn.zeros(4,4)

  iter <- 1
  while (iter <= con$itermax & bestval > con$VTR) {     
    popold <- pop ## save old population
    
    ind <- sample(1:4) ## index pointer array

    a1 <- sample(1:NP) ## shuffle locations and rotate vectors
    rt <- (rot + ind[1]) %% NP 
    a2 <- a1[rt + 1] 
    rt <- (rot + ind[2]) %% NP 
    a3 <- a2[rt + 1]
    rt <- (rot + ind[3]) %% NP
    a4 <- a3[rt + 1]     
    rt <- (rot + ind[4]) %% NP
    a5 <- a4[rt + 1]
      
    pm1 <- popold[a1,] ## shuffled populations 1 - 5
    pm2 <- popold[a2,]
    pm3 <- popold[a3,] 
    pm4 <- popold[a4,] 
    pm5 <- popold[a5,] 
    
    bm <- matrix(rep(bestmemit, NP), nrow = NP, byrow = TRUE) ## population filled with
    ## the best member of the last iteration
      
    mui <- matrix(runif(NP * d), nrow = NP) < con$CR ## all random numbers < CR are 1, 0 otherwise
      
    if (con$strategy > 5)
      st <- con$strategy - 5 ## binomial crossover
    else {
      st <- con$strategy ## exponential crossover
      mui <- sort(t(mui)) ## transpose, collect 1's in each column
        
      for (i in 1:NP) {
        n <- floor(runif(1) * d)
        if (n > 0) {
          rtd <- (rotd + n) %% d
          mui[,i] <- mui[rtd + 1,i] ## rotate column i by n 
        }
      }
      mui <- t(mui) ## transpose back
    }
    
    mpo <- mui < 0.5 

    if (st == 1) { ## best / 1
      ui <- bm + con$F * (pm1 - pm2) ## differential variation
      ui <- popold * mpo + ui * mui ## crossover
    }
    else if (st == 2) { ## rand / 1
      ui <- pm3 + con$F * (pm1 - pm2) ## differential variation
      ui <- popold * mpo + ui * mui ## crossover
    }
    else if (st == 3) { ## rand-to-best / 1
      ui <- popold + con$F * (bm - popold) + con$F * (pm1 - pm2) ## differential variation
      ui <- popold * mpo + ui * mui ## crossover
    }
    else if (st == 4) { ## best / 2
      ui <- bm + con$F * (pm1 - pm2 + pm3 - pm4) ## differential variation
      ui <- popold * mpo + ui * mui ## crossover
    }
    else { ## rand / 2                
      ui <- pm5 + con$F * (pm1 - pm2 + pm3 - pm4) ## differential variation
      ui <- popold * mpo + ui * mui ## crossover
    }

    for (i in 1:NP)
      ui[i,] <- fn.checkBoundaries(ui[i,], lower, upper) ## check whether
    ## the components are within the boundaries

    nfeval <- nfeval + NP
    tempval <- apply(ui, 1, FUN, ...) ## check cost of competitor
    ichange <- tempval <= val
    val[ichange] <- tempval[ichange]
    pop[ichange,] <- ui[ichange,]
    bestval <- bestvalit <- min(val)
    ibest <- match(bestvalit, val)
    bestmem <- bestmemit <- pop[ibest,]
    
    if (con$refresh > 0 & (iter %% con$refresh) == 0) { 
      cat("iteration: ", iter,
          "bestmem: " , bestmem,
          "bestval: ", bestval, "\n")
    }
    iter <- iter + 1
  }

  ##
  ## Outputs
  ##
  
  r <- list(bestmem, bestval, nfeval)
  names(r) <- c("bestmem", "bestval", "nfeval")
  return(r)
}

