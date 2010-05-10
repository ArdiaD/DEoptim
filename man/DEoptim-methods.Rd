\name{DEoptim-methods}
\alias{DEoptim-methods}
\alias{plot.DEoptim}
\alias{summary.DEoptim}
\title{DEoptim-methods}
\description{Methods for DEoptim objects.}
\usage{
\method{summary}{DEoptim}(object, \dots)
\method{plot}{DEoptim}(x, plot.type = c("bestmemit", "bestvalit", "storepop"), \dots)
}
\arguments{
  \item{object}{An object of class \code{DEoptim}; usually, a result
    of a call to \code{\link{DEoptim}}.}
  \item{x}{An object of class \code{DEoptim}; usually, a result
    of a call to \code{\link{DEoptim}}.}
  \item{plot.type}{Should we plot the best member at each iteration, the best value
    at each iteration or the intermediate populations?}
  \item{\dots}{Further arguments passed to or from other methods.}
}
\note{
  Please cite the package in publications. Use \code{citation("DEoptim")}.
}
\author{
  David Ardia \email{david.ardia@unifr.ch} and 
  Katharine Mullen \email{katharine.mullen@nist.gov}.
}
\seealso{
  \code{\link{DEoptim}} and \code{\link{DEoptim-methods}}.
}
\examples{
  ## Rosenbrock Banana function
  Rosenbrock <- function(x){
    x1 <- x[1]
    x2 <- x[2]
    100 * (x2 - x1 * x1)^2 + (1 - x1)^2
  }

  lower <- c(-10, -10)
  upper <- -lower
  
  set.seed(1234)

  outDEoptim <- DEoptim(Rosenbrock, lower, upper)
  
  ## print output information
  summary(outDEoptim)

  ## plot the best members
  plot(outDEoptim, type = 'b')

  ## plot the best values
  dev.new()
  plot(outDEoptim, plot.type = "bestvalit", type = 'b', col = 'blue')

  ## rerun the optimization, and store intermediate populations
  outDEoptim <- DEoptim(Rosenbrock, lower, upper,
                        DEoptim.control(itermax = 500,
                        storepopfrom = 1, storepopfreq = 2))
  summary(outDEoptim)
  
  ## plot intermediate populations
  dev.new()
  plot(outDEoptim, plot.type = "storepop")

  ## Wild function
  Wild <- function(x)
    10 * sin(0.3 * x) * sin(1.3 * x^2) +
       0.00001 * x^4 + 0.2 * x + 80

  outDEoptim = DEoptim(Wild, lower = -50, upper = 50,
                       DEoptim.control(trace = FALSE, storepopfrom = 50,
                       storepopfreq = 1))
  
  plot(outDEoptim, type = 'b')

  dev.new()
  plot(outDEoptim, plot.type = "bestvalit", type = 'b')

\dontrun{
  ## an example with a normal mixture model: requires package mvtnorm
  library(mvtnorm)

  ## neg value of the density function
  negPdfMix <- function(x) {
     tmp <- 0.5 * dmvnorm(x, c(-3, -3)) + 0.5 * dmvnorm(x, c(3, 3))
     -tmp
  }

  ## wrapper plotting function
  plotNegPdfMix <- function(x1, x2)
     negPdfMix(cbind(x1, x2))

  ## contour plot of the mixture
  x1 <- x2 <- seq(from = -10.0, to = 10.0, by = 0.1)
  thexlim <- theylim <- range(x1)
  z <- outer(x1, x2, FUN = plotNegPdfMix)
  
  contour(x1, x2, z, nlevel = 20, las = 1, col = rainbow(20),
     xlim = thexlim, ylim = theylim)

  set.seed(1234)

  outDEoptim <- DEoptim(negPdfMix, c(-10, -10), c(10, 10),
     DEoptim.control(NP = 100, itermax = 100, storepopfrom = 1,
     storepopfreq = 5))

  ## convergence plot
  dev.new()
  plot(outDEoptim)
  
  ## the intermediate populations indicate the bi-modality of the function
  dev.new()
  plot(outDEoptim, plot.type = "storepop")
}
}
\keyword{methods}
