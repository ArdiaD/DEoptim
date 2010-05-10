\name{DEoptim-methods}
\alias{DEoptim-methods}
\alias{plot.DEoptim}
\alias{summary.DEoptim}
\title{DEoptim-methods}
\description{
  Methods for DEoptim objects.
}
\usage{
\method{summary}{DEoptim}(object, \dots)
\method{plot}{DEoptim}(x, plot.type = c("bestmemit","bestvalit","storepop"), \dots)
}
\arguments{
  \item{object, x}{An object of class~\code{DEoptim}; usually, a result
    of a call to~\code{\link{DEoptim}}.}
  \item{plot.type}{Should we plot the best member at each iteration, the best value 
    at each iteration or the intermediate populations?}
  \item{\dots}{Further arguments passed to or from other methods.}
}
\note{
  Please cite the package in publications. Use \code{citation("DEoptim")}.
}
\author{
  David Ardia~\email{david.ardia@unifr.ch}
}
\examples{
  ## Rosenbrock Banana function
  Rosenbrock <- function(x){
    x1 <- x[1]
    x2 <- x[2]
    100 * (x2 - x1 * x1)^2 + (1 - x1)^2
  }

  lower <- c(-10,-10)
  upper <- -lower
  r <- DEoptim(Rosenbrock, lower, upper)
  summary(r) ## print summary of the DEoptim object

  par(mfrow = c(1,2))
  plot(r, type = 'b') ## plot the best members
  plot(r, plot.type = "bestvalit", type = 'b', col = 'blue') ## plot the best values

  ## rerun the optimization, and store intermediate populations
  r <- DEoptim(Rosenbrock, lower, upper, control=list(NP=400, itermax=30, storepopfrom=1, storepopfreq=2))
  summary(r)
  
  par(mfrow=c(1,1))
  plot(r, plot.type = "storepop") ## plot intermediate populations
}
\keyword{methods}
