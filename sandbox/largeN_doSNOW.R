#############################################################################
## Set up objective function to minimize 
#############################################################################

library(PerformanceAnalytics)
library(DEoptim)
library(doSNOW)
library(parallel)

load("10y_returns.rda")
load("random_portfolios.rda")

rng <- 1:30
R <- R[,rng]
rp <- rp[1:300,rng]

m3 <- PerformanceAnalytics:::M3.MM(R)
m4 <- PerformanceAnalytics:::M4.MM(R)

mu <- colMeans(R)
sigma <- cov(R)
N <- ncol(R)
lower <- rep(0,N) 
upper <- rep(1,N)

# objective function (do not specify non-varying parameters in function definition)
obj <- function(w) {
  if(sum(w)==0) w <- w + 1e-2
  w <- w / sum(w)
  CVaR <- ES(weights=w, method="modified",
    portfolio_method="component", mu=mu, sigma=sigma, m3=m3, m4=m4)
  tmp1 <- CVaR$MES
  tmp2 <- max(CVaR$pct_contrib_MES - 0.05, 0)
  out <- tmp1 + 1e3 * tmp2
  return(out)
}

#############################################################################
## Call DEoptim on one CPU 
#############################################################################

## Note that parallelType = 0 means that only one CPU will be utilized,
## which was the default behavior of DEoptim prior to version 2-2.0 

controlDE <- list(NP = nrow(rp), initialpop = rp, trace = 1, itermax = 5,
                  reltol = 0.000001, steptol = 150, c = 0.4, strategy = 6,
                  parallelType = 0)

## Call DEoptim using only one CPU
set.seed(1234)
timeONECORE <- system.time(out1 <- DEoptim(fn=obj, lower=lower, upper=upper,
  control=controlDE))
out1$optim$iter
out1$optim$bestval

#############################################################################
## Call DEoptim on multiple CPU's, using the parallel package   
#############################################################################

## Note that the packages needed are listed in packages
## Note that packVar gives the arguments needed 

controlDE1 <- list(NP = nrow(rp), initialpop = rp, trace = 1, itermax = 5,
                   reltol = 0.000001, steptol = 150, c = 0.4, strategy = 6,
                   parallelType = 1, 
                   packages = list("PerformanceAnalytics"),
                   parVar = list("mu", "sigma", "m3", "m4"))

set.seed(1234)

timeALLCORES1 <- system.time(out2 <- DEoptim(fn=obj, lower=lower, upper=upper,
  control=controlDE1))

############################################################################
## Call DEoptim on multiple CPU's, using the foreach package 
#############################################################################

## Find out how many CPUs are available on the local machine 
nC <- detectCores() 

## Get ready to evaluate on all these cores 
cl <- makeSOCKcluster(nC)

## load any necessary packages 
clusterEvalQ(cl, library(PerformanceAnalytics)) 

## copy any necessary objects
clusterExport(cl, list("mu", "sigma", "m3", "m4"))

## register foreach backend
registerDoSNOW(cl) 

controlDE <- list(
                  NP = nrow(rp), initialpop = rp, trace = 1, itermax = 5,
                  reltol = 0.000001, steptol = 150, c = 0.4, strategy = 6,
                  parallelType = 2)

## Call DEoptim using all available cores 
## *do not* pass mu, sigma, m3, m4 
set.seed(1234)
timeALLCORES2 <- system.time(out2 <- DEoptim(fn=obj, lower=lower, upper=upper,
  control=controlDE)) 
## stop cluster
stopCluster(cl) 

## Compare timings
timeONECORE
timeALLCORES1
timeALLCORES2


