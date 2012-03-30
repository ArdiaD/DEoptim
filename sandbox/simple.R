#############################################################################
## Set up objective function to minimize 
#############################################################################

library(DEoptim) 

Genrose <- function(x) {        ## One generalization of the Rosenbrock banana valley function (n parameters)
  n <- length(x)
  ## make it take some time ... 
  Sys.sleep(.001) 
  1.0 + sum (100 * (x[-n]^2 - x[-1])^2 + (x[-1] - 1)^2)
}

# get some run-time on simple problems
maxIt <- 250                     
n <- 5

oneCore <- system.time( DEoptim(fn=Genrose, lower=rep(-25, n), upper=rep(25, n),
                   control=list(NP=10*n, itermax=maxIt)))

withParallel <-  system.time( DEoptim(fn=Genrose, lower=rep(-25, n), upper=rep(25, n),
                   control=list(NP=10*n, itermax=maxIt, parallelType=1)))

oneCore
withParallel
