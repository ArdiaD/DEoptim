library(DEoptim)

n <- 5
lims <- rep(25, n)
ctrl <- list(NP = 10*n, itermax = 500, trace = 0)

## One generalization of the Rosenbrock banana valley function (n parameters)
genrose <- function(x, a, b) {
  n <- length(x)
  1.0 + sum(b * (x[-n]^2 - x[-1])^2 + (x[-1] - a)^2)
}
Rf <- function() {
  DEoptim(genrose, -lims, lims, ctrl, a = 1.0, b = 10.0)
}

fn_ptr <- getNativeSymbolInfo("genrose", "DEoptim")
Cf <- function() {
  DEoptim(fn_ptr, -lims, lims, ctrl, a = 1.0, b = 10.0)
}

system.time(Rf())
system.time(Cf())

library(microbenchmark)
microbenchmark(Cf(), Rf())

set.seed(21); outCf <- Cf();
set.seed(21); outRf <- Rf();
all.equal(outCf, outRf)
