
library(PerformanceAnalytics)
library(PortfolioAnalytics)
library(DEoptim)
library(doSNOW)

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

controlDE <- list(
  NP = nrow(rp), initialpop = rp, trace = 1, itermax = 5,
  reltol = 0.000001, steptol = 150, c = 0.4, strategy = 6 )

set.seed(1234)
system.time(out1 <- DEoptim(fn=obj, lower=lower, upper=upper,
  control=controlDE))
out1$optim$iter
out1$optim$bestval


cl <- makeSOCKcluster(2)
clusterEvalQ(cl, library(PerformanceAnalytics)) # load any necessary libraries
clusterExport(cl, list("mu", "sigma", "m3", "m4")) # copy any necessary objects
registerDoSNOW(cl) # register foreach backend
set.seed(1234)
system.time(out2 <- DEoptim(fn=obj, lower=lower, upper=upper,
  control=controlDE)) # *do not* pass mu, sigma, m3, m4
stopCluster(cl) # stop cluster
out2$optim$iter
out2$optim$bestval

