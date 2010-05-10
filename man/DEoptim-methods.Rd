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
\method{plot}{DEoptim}(x, plot.type = c("bestmemit","bestvalit"), \dots)
}
\arguments{
  \item{object, x}{An object of class~\code{DEoptim}; usually, a result
    of a call to~\code{\link{DEoptim}}.}
  \item{plot.type}{Should we plot the best member at each iteration or
    the best value at each iteration?}
  \item{\dots}{Further arguments passed to or from other methods.}
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
  plot(r, plot.type = "bestvalit", type = 'b', col = 'blue') ## and the best values
}
\keyword{methods}
