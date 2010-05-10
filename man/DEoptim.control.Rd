\name{DEoptim.control}
\alias{DEoptim.control}
\title{Control various aspects of the DEoptim implementation}
\description{
  Allow the user to set some characteristics of the 
  Differential Evolution optimization algorithm implemented
  in \code{DEoptim}.
}
\usage{
DEoptim.control(VTR = -Inf, strategy = 2, bs = TRUE, NP = 50,
   itermax = 200, CR = 0.5, F = 0.8, trace = TRUE,
   initialpop = NULL, storepopfrom = itermax + 1,
   storepopfreq = 1)
}
\arguments{
   \item{VTR}{The value to be reached. The optimization process
      will stop if either the maximum number of iterations \code{itermax}
      is reached or the best parameter vector \code{bestmem} has found a value
      \code{fn(bestmem) <= VTR}. Default to \code{-Inf}.}
  \item{strategy}{Defines the Differential Evolution
    strategy used in the optimization procedure:\cr
     \code{1}: DE / rand / 1 / bin (classical strategy)\cr
     \code{2}: DE / local-to-best / 1 / bin \cr
     \code{3}: DE / best / 1 / bin with jitter\cr
     \code{4}: DE / rand / 1 / bin with per-vector-dither\cr
     \code{5}: DE / rand / 1 / bin with per-generation-dither\cr
     any value not above:
     variation to DE / rand / 1 / bin: either-or-algorithm.
    }
    \item{bs}{Enables best of parent and child selection if \code{TRUE} or
      DEoptim standard trial vs. target selection if \code{FALSE}.
      Default is \code{TRUE}.}
    \item{NP}{Number of population members. Default to \code{50}.
      Cannot be set larger than 2000.}
    \item{itermax}{The maximum iteration (population generation) allowed.
      Default is \code{200}.}
    \item{CR}{Crossover probability from interval [0,1]. Default
      to \code{0.5}.}
    \item{F}{Stepsize from interval [0,2]. Default to \code{0.8}.}
    \item{trace}{Printing of progress occurs? Default to \code{TRUE}.}
    \item{initialpop}{An initial population used as a starting
      population in the optimization procedure. Maybe useful to speed up
      the convergence. Default to \code{NULL}.}
    \item{storepopfrom}{From which population should the following
      intermediate populations be stored in memory. Default to
      \code{itermax+1}, i.e., no intermediate population is stored.}
    \item{storepopfreq}{The frequency of populations'
      storage. Default to \code{1}, i.e. every intermediate population
      is memorized.}
  }
  \value{
    A \code{list} with components:
    \item{VTR}{}
    \item{strategy}{}
    \item{bs}{}
    \item{NP}{}
    \item{itermax}{}
    \item{CR}{}
    \item{F}{}
    \item{trace}{}
    with meanings as explained under \sQuote{Arguments}.
}
\references{
  Differential Evolution homepage:
  \url{http://www.icsi.berkeley.edu/~storn/code.html}
  
  Price, K.V., Storn, R.M., Lampinen J.A. (2005).
  \emph{Differential Evolution - A Practical Approach to Global Optimization}.
  Springer-Verlag. ISBN 3540209506.
}
\seealso{
  \code{\link{DEoptim}} and \code{\link{DEoptim-methods}}.
}
\note{
  Please cite the package in publications. Use \code{citation("DEoptim")}.
}
\author{
  David Ardia \email{david.ardia@unifr.ch} and 
  Katharine Mullen \email{katharine.mullen@nist.gov}.
}
\examples{
DEoptim.control(NP = 20)

DEoptim.control(NP = 20, itermax = 100, trace = FALSE)
}
\keyword{nonlinear}
\keyword{optimize}

