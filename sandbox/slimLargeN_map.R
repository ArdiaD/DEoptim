setwd("~/R/packages/deoptim/pkg/DEoptim/sandbox")
#suppressMessages({
library(PerformanceAnalytics)
library(DEoptim)
#})

load("10y_returns.rda")
load("random_portfolios.rda")

rng <- 1:20
R <- R[,rng]
rp <- rp[1:200,rng]

mu <- colMeans(R)
sigma <- cov(R)
N <- ncol(R)
lower <- rep(0,N) 
upper <- rep(1,N)

obj <- function(w) {
  if(sum(w)==0) w <- w + 1e-2
  w <- w / sum(w)
  CVaR <- ES(weights=w, method="gaussian",
    portfolio_method="component", mu=mu, sigma=sigma)
  tmp1 <- CVaR$ES
  tmp2 <- max(CVaR$pct_contrib_ES - 0.05, 0)
  out <- tmp1 + 1e3 * tmp2
  return(out)
}

controlDE <- list(
  NP = nrow(rp), initialpop = rp, trace = 10, itermax = 50,
  reltol = 0.000001, steptol = 150, c = 0.4, strategy = 6 )

set.seed(1234) # reset random chain before DE step
system.time(out1 <- DEoptim(fn=obj, lower=lower, upper=upper, control=controlDE))
out1$optim$iter
out1$optim$bestval

# this should throw a ton of warnings
set.seed(1234)
system.time(out <- DEoptim(fn=obj, lower=lower, upper=upper,
  control=controlDE, fnMap=function(x) 1:length(x)))
out$optim$iter
out$optim$bestval

# trivial full-investment function
mappingFun <- function(x) {
  x[which(order(x) < 6)] <- 0
  x <- round(x,2) # produce some dups
  x/sum(x)
}

set.seed(1234)
system.time(out <- DEoptim(fn=obj, lower=lower,
  upper=upper, control=controlDE, fnMap=mappingFun))
out$optim$iter
out$optim$bestval
