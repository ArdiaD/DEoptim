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
    return a scalar result.}
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
      \code{FUN(bestmem) <= VTR}. Defaults to \code{10e-6}.}
    \item{\code{itermax}}{The maximum iteration (population generation) allowed.
      Defaults is \code{200}.}
    \item{\code{NP}}{Number of population members. Defaults to \code{50}.}
    \item{\code{F}}{Stepsize from interval [0,2]. Defaults to \code{0.8}.}
    \item{\code{CR}}{Crossover probability from interval [0,1]. Defaults
      to \code{0.5}.}
    \item{\code{strategy}}{Defines the DE-strategy used in the
      optimization process:
      \describe{
        \item{\code{1}}{best/1/exponential}
        \item{\code{2}}{rand/1/exponential}
        \item{\code{3}}{rand-to-best/1/exponential}
        \item{\code{4}}{best/2/exponential}
        \item{\code{5}}{rand/2/exponential}
        \item{\code{6}}{best/1/binomial}
        \item{\code{7}}{rand/1/binomial}
        \item{\code{8}}{rand-to-best/1/binomial}
        \item{\code{9}}{best/2/binomial}
        \item{\code{10}}{rand/2/binomial}
      }
      Experiments suggest that binomial likes to have a slightly larger \code{CR}
      than exponential. By default \code{strategy} is \code{7}.}
    \item{\code{refresh}}{The frequency of reports. Defaults to
      every \code{10} iterations.}
  }
}
\value{
  A list with components:
    \item{bestmem}{The best set of parameters found.}
    \item{bestval}{The value of \code{FUN} corresponding to \code{bestmem}.}
    \item{nfeval}{Number of function evaluations.}
}
\note{
  \code{DEoptim} is a R-vectorized variant of the Differential Evolution algorithm
  initialy developed by R. Storn \email{storn@icsi.berkeley.edu},
  International Computer Science Institute (ICSI), 1947 Center Street, Suite 600,
  Berkeley, CA 94704. \url{http://http.icsi.berkeley.edu/~storn}.

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
  
  If you choose binomial crossover like, rand/1/binomial, \code{CR} is usually 
  higher than in the exponential crossover variant (in this example 
  rand/1/exponential). Still, \code{F} and \code{CR} are both generally in 
  the range [0.5,1] for most problems we have encountered. Keep in mind that 
  different problems usually require different settings for \code{NP}, \code{F} 
  and \code{CR}. 
  
  If you still get misconvergence you might want to try a different method. 
  We mostly use rand/1/... or best/2/... (for the latter \code{F=0.5} is the 
  standard choice). The crossover method is not so important although Ken Price 
  claims that binomial is never worse than exponential. In case of misconvergence 
  also check your choice of objective function. There might be a better one to 
  describe your problem. Any knowledge that you have about the problem should be 
  worked into the objective function. A good objective function can make all 
  the difference.

  The R-adaptation \code{DEoptim} has properties which differ from the 
  original DE version:
  \describe{
    \item{1.}{The random selection of vectors is performed by shuffling the
      population array. Hence a certain vector cannot be chosen twice in the same
      term of the perturbation expression.}
    \item{2.}{Due to the vectorized expressions \code{DEoptim} executes fairly fast.}
    \item{3.}{The parameters are constrained within boundaries.}
  }
}

\references{
  \url{http://www.icsi.berkeley.edu/~storn}

  \url{http://www2.lut.fi/~jlampine/debiblio.htm}
  
  Nocedal, J. and Wright, S. J. (1999) \emph{Numerical Optimization}.
  Springer.

  Price, K., Storn, R. and J. Lampingen (2005) \emph{Differential Evolution
  - A Practical Approach to Global Optimization over Continuous Spaces}.
  Springer.

  Storn, R. and Price, K. (1995) Differential Evolution - a Simple and Efficient
  Adaptive Scheme for Global Optimization over Continuous Spaces.
  \emph{Technical Report TR-95-012}.
}

\author{David Ardia}

\seealso{\code{\link{nlm}}, \code{\link{optim}}, \code{\link{constrOptim}}}

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
}
\keyword{nonlinear}
\keyword{optimize}
