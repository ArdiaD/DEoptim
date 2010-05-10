\name{DEoptim}
\alias{DEoptim}
\title{Differential Evolution Optimization}
\concept{minimization}
\description{
  Performs evolutionary optimization via the Differential Evolution algorithm.
}
\usage{
DEoptim(FUN, lower, upper, control = list(), ...)
}
\arguments{
  \item{FUN}{A function to be minimized, with first argument the vector
    of parameters over which minimization is to take place. It should
    return a scalar result. \code{NA} and \code{NaN} values are not allowed.}
  \item{lower, upper}{Bounds on the variables.}
  \item{control}{A list of control parameters. See *Details*.}
  \item{...}{Further arguments to be passed to \code{FUN}.}
}
\details{
  \code{DEoptim} performs minimization of \code{FUN}.

  The \code{control} argument is a list that can supply any of
  the following components:

  \describe{
    \item{\code{VTR}}{The value to be reached. The optimization process
      will stop if either the maximum number of iterations \code{itermax}
      is reached or the best parameter vector \code{bestmem} has found a value
      \code{FUN(bestmem) <= VTR}. Default to \code{-Inf}.}
    \item{\code{itermax}}{The maximum iteration (population generation) allowed.
      Default is \code{200}.}
    \item{\code{NP}}{Number of population members. Default to \code{50}.}
    \item{\code{F}}{Stepsize from interval [0,2]. Default to \code{0.8}.}
    \item{\code{CR}}{Crossover probability from interval [0,1]. Default
      to \code{0.5}.}
    \item{\code{initial}}{An initial population used as a starting
      population in the optimization procedure. Maybe useful to speed up
      the convergence. Defaults to \code{NULL}.}
    \item{\code{strategy}}{Defines the binomial DE-strategy used in the
      optimization procedure:
      \describe{
        \item{\code{1}}{best/1}
        \item{\code{2}}{rand/1}
        \item{\code{3}}{rand-to-best/1}
        \item{\code{4}}{best/2}
        \item{\code{5}}{rand/2}
      }
      
      By default \code{strategy} is \code{2}. See references below for details.}

    \item{\code{refresh}}{The frequency of reports. Default to
      every \code{10} iterations.}
    \item{\code{digits}}{The number of digits to print when printing numeric values at each
      iteration. Default to 4.}
  }
}
\value{
  A list of lists of the class \code{DEoptim}.\cr

  list \code{optim} contains the followings:\cr
  \code{bestmem}: the best set of parameters found.\cr
  \code{bestval}: the value of \code{FUN} corresponding to \code{bestmem}.\cr
  \code{nfeval}: number of function evaluations.\cr
  \code{iter}: number of procedure iterations.\cr

  list \code{member} contains the followings:\cr
  \code{lower}: the lower boundary.\cr
  \code{upper}: the upper boundary.\cr
  \code{bestvalit}: the best value of \code{FUN} at each iteration.\cr
  \code{bestmemit}: the best member at each iteration.\cr
  \code{pop}: the population generated at the last iteration.\cr
}
\note{
  \code{DEoptim} is a \R-vectorized variant of the Differential Evolution algorithm
  initialy developed by Rainer Storn \email{storn@icsi.berkeley.edu},
  International Computer Science Institute (ICSI), 1947 Center Street, Suite 600,
  Berkeley, CA 94704.
  
  If you experience misconvergence in the optimization process you usually
  have to increase the value for \code{NP}, but often you only have to adjust
  \code{F} to be a little lower or higher than \code{0.8}. If you increase
  \code{NP} and simultaneously lower \code{F} a little, convergence is more
  likely to occur but generally takes longer, i.e. \code{DEoptim} is getting
  more robust (there is always a convergence speed/robustness tradeoff).

  \code{DEoptim} is much more sensitive to the choice of \code{F} than it is to
  the choice of \code{CR}. \code{CR} is more like a fine tuning element. High
  values of \code{CR} like \code{CR=1} give faster convergence if convergence
  occurs. Sometimes, however, you have to go down as much as \code{CR=0} to
  make \code{DEoptim} robust enough for a particular problem.

  The \R-adaptation \code{DEoptim} has properties which differ from the
  original DE version:
  \describe{
    \item{1.}{The random selection of vectors is performed by shuffling the
      population array. Hence a certain vector cannot be chosen twice in the same
      term of the perturbation expression.}
    \item{2.}{Due to the vectorized expressions \code{DEoptim} executes fairly fast.}
    \item{3.}{The parameters are constrained within boundaries.}
    \item{4.}{An initial population can be given as a starting point for
      the DE-optimization. This may speed up the convergence if the
      optimization procedure has to be run many times for sligthly
      different data sets.}
  }
  
  To perform a maximization (instead of minimization) of a given
  function, simply define a new function which is the opposite of the
  function to maximize and apply \code{DEoptim} to it.
  
  To integrate additional constraints on the parameters \code{x} of
  \code{FUN(x)}, for instance \code{x[1] + x[2]^2 < 2}, integrate the
  constraint within the function to optimize, for instance: 
  \preformatted{
    FUN <- function(x)\{
      if (x[1] + x[2]^2 < 2)\{
	r <- Inf
      else\{
	...
      \}
      return(r)
    \}
  }

  Note that \code{DEoptim} stops if any \code{NA} or \code{NaN} value is
  obtained. You have to redefine your function to handle these values
  (for instance, set \code{NA} to \code{Inf} in your objective function).
}
\references{
  Differential Evolution homepage :

  \url{http://www.icsi.berkeley.edu/~storn/code.html}

  Some useful books:
  
  Price, K. V., Storn, R. M. and Lampinen J. A. (2005) \emph{Differential Evolution
    - A Practical Approach to Global Optimization}.
  Springer. ISBN 3540209506.

  Nocedal, J. and Wright, S. J. (1999) \emph{Numerical Optimization}.
  Springer. ISBN 0387987932.
}

\author{David Ardia \email{david.ardia@unifr.ch} for the \R-port;
  Rainer Storn \email{storn@icsi.berkeley.edu} for the Differential Evolution algorithm.}

\seealso{\code{\link{DEoptim-methods}} for methods on \code{DEoptim} object;
  \code{\link{optim}} or \code{\link{constrOptim}} for constrained optimization.}

\examples{
  ## Rosenbrock Banana function
  Rosenbrock <- function(x){
    x1 <- x[1]
    x2 <- x[2]
    100 * (x2 - x1 * x1)^2 + (1 - x1)^2
  }

  lower <- c(-10,-10)
  upper <- -lower
  DEoptim(Rosenbrock, lower, upper)
  DEoptim(Rosenbrock, lower, upper, 
    control = list(NP = 100, refresh = 1))
  DEoptim(Rosenbrock, lower, upper, 
    control = list(NP = 50, itermax = 200, F = 1.5, 
    CR = 0.2, refresh = 1))
  DEoptim(Rosenbrock, lower, upper,
    control = list(NP = 80, itermax = 400, F = 1.2,
    CR = 0.7, refresh = 1))

  ## 'Wild' function, global minimum at about -15.81515
  Wild <- function(x)
    10 * sin(0.3*x) * sin(1.3*x^2) + 
       0.00001 * x^4 + 0.2 * x + 80
  plot(Wild, -50, 50, n = 1000, 
    main = "DEoptim minimizing 'Wild function'")
  DEoptim(Wild, lower = -50, upper = 50,
      control = list(NP = 50, refresh = 1))
  DEoptim(Wild, lower = -50, upper = 50, 
      control = list(NP = 50, refresh = 1, digits = 8))
}
\keyword{nonlinear}
\keyword{optimize}
