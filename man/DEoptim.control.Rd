\name{DEoptim.control}
\alias{DEoptim.control}
\title{Control various aspects of the DEoptim implementation}
\description{
  Allow the user to set some characteristics of the 
  differential evolution optimization algorithm implemented
  in \code{DEoptim}.
}
\usage{
DEoptim.control(VTR = -Inf, strategy = 2, bs = FALSE, NP = 50,
   itermax = 200, CR = 0.5, F = 0.8, trace = TRUE,
   initialpop = NULL, storepopfrom = itermax + 1,
   storepopfreq = 1, checkWinner = FALSE, avWinner = TRUE)
}
\arguments{
   \item{VTR}{The value to be reached. The optimization process
      will stop if either the maximum number of iterations \code{itermax}
      is reached or the best parameter vector \code{bestmem} has found a value
      \code{fn(bestmem) <= VTR}. Default to \code{-Inf}.}
  \item{strategy}{Defines the differential evolution
    strategy used in the optimization procedure:\cr
     \code{1}: DE / rand / 1 / bin (classical strategy)\cr
     \code{2}: DE / local-to-best / 1 / bin \cr
     \code{3}: DE / best / 1 / bin with jitter\cr
     \code{4}: DE / rand / 1 / bin with per-vector-dither\cr
     \code{5}: DE / rand / 1 / bin with per-generation-dither\cr
     any value not above:
     variation to DE / rand / 1 / bin: either-or-algorithm.  Default
     strategy is currently \code{2}.
    }
    \item{bs}{If \code{FALSE} then every mutant will be tested against a
      member in the previous generation, and the best value will proceed
      into the next generation (this is standard trial vs. target
      selection).   If \code{TRUE} then the old generation and \code{NP}
      mutants will be sorted by their associated objective function
      values, and the best \code{NP} vectors will proceed into the next
      generation (best of parent and child selection). Default is
    \code{FALSE}.}
    \item{NP}{Number of population members. Defaults to \code{50}.  For
      many problems it is best to set
      \code{NP} to be at least 10 times the length
      of the parameter vector.  }
    \item{itermax}{The maximum iteration (population generation) allowed.
      Default is \code{200}.}
    \item{CR}{Crossover probability from interval [0,1]. Default
      to \code{0.9}.}
    \item{F}{Stepsize from interval [0,2]. Default to \code{0.8}.}
    \item{trace}{Printing of progress occurs? Default to \code{TRUE}.}
    \item{initialpop}{An initial population used as a starting
      population in the optimization procedure. May be useful to speed up
      the convergence. Default to \code{NULL}.  If given, each member of
    the initial population should be given as a row of a numeric matrix, so that
   \code{initialpop} is a matrix with \code{NP} rows and a number of
   columns equal to the length of the parameter vector to be optimized. }
    \item{storepopfrom}{From which generation should the following
      intermediate populations be stored in memory. Default to
      \code{itermax+1}, i.e., no intermediate population is stored.}
    \item{storepopfreq}{The frequency with which populations are stored.
      Default to \code{1}, i.e. every intermediate population
      is stored.}
    \item{checkWinner}{Logical value indicating whether to re-evaluate
      the objective function using the winning parameter vector if this
      vector remains the same between
    generations.  This may be useful for the optimization of a noisy
    objective function.  If \code{checkWinner=TRUE} and \code{avWinner=FALSE}
   then the  value associated with re-evaluation of
   the objective function is used in the next generation.
   Default to \code{FALSE}. }
 \item{avWinner}{Logical value.  If  \code{checkWinner=TRUE} and
   \code{avWinner=TRUE} then the objective
   function value associated with the winning member represents the
   average of all evaluations of the objective function over the course
   of the `winning streak' of the best population member. This option
   may be useful for optimization of noisy objective functions, and is
   interpretted only if \code{checkWinner=TRUE}.
   The default value is \code{TRUE}.}
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
    \item{initialpop}{}
    \item{storepopfrom}{}
    \item{storepopfreq}{}
    \item{checkWinner}{}
    \item{avWinner}{}
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

