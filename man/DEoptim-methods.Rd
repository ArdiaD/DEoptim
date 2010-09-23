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
  \item{object}{an object of class \code{DEoptim}; usually, a result
    of a call to \code{\link{DEoptim}}.}
  \item{x}{an object of class \code{DEoptim}; usually, a result
    of a call to \code{\link{DEoptim}}.}
  \item{plot.type}{should we plot the best member at each iteration, the best value
    at each iteration or the intermediate populations?}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
  Members of the class \code{DEoptim} have a \code{plot} method that
  accepts the argument \code{plot.type}. \code{plot.type = "bestmemit"} results
  in a plot of the parameter values that represent the lowest value of the objective function
  each generation. \code{plot.type = "bestvalit"} plots the best value of
  the objective function each generation. Finally, \code{plot.type = "storepop"} results in a plot of
  stored populations (which are only available if these have been saved by
  setting the \code{control} argument of \code{DEoptim} appropriately). Storing intermediate populations 
  allows us to examine the progress of the optimization in detail.   
  A summary method also exists and returns the best parameter vector, the best value of the objective function,
  the number of generations optimization ran, and the number of times the 
  objective function was evaluated. 
}
\note{
  Further details and examples of the \R package \pkg{DEoptim} can be found
  in Mullen et al. (2009) and Ardia et al. (2010).

  Please cite the package in publications. Use \code{citation("DEoptim")}.
}
\author{
  David Ardia, Katharine Mullen \email{katharine.mullen@nist.gov}, 
  Brian Peterson and Joshua Ulrich.
}
\seealso{
  \code{\link{DEoptim}} and \code{\link{DEoptim.control}}.
}
\references{
  Mullen, K.M., Ardia, D., Gil, D.L, Windover, D., Cline, J. (2009)
  \pkg{DEoptim}: An \R Package for Global Optimization by Differential Evolution.
  URL \url{http://ssrn.com/abstract=1526466}

  Ardia, D., Boudt, K., Carl, P., Mullen, K.M., Peterson, B.G. (2010)
  Differential Evolution (\pkg{DEoptim}) for Non-Convex Portfolio Optimization.
  URL \url{http://ssrn.com/abstract=1584905}
}
\examples{
  ## Rosenbrock Banana function
  ## The function has a global minimum f(x) = 0 at the point (0,0).  
  ## Note that the vector of parameters to be optimized must be the first 
  ## argument of the objective function passed to DEoptim.
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
