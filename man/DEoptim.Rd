\name{DEoptim}
\alias{DEoptim}
\title{Differential Evolution Optimization}
\concept{minimization}
\description{
  Performs evolutionary global optimization via the Differential Evolution algorithm.
}
\usage{
DEoptim(fn, lower, upper, control = DEoptim.control(), ...)
}
\arguments{
  \item{fn}{the function to be optimized (minimized). The function should have as its first 
    argument the vector of real-valued parameters to optimize, and return a scalar real result. \code{NA} 
    and \code{NaN} values are not allowed.}
  \item{lower, upper}{two vectors specifying scalar real lower and upper bounds on each parameter to be optimized, so that the i-th element 
    of \code{lower} and \code{upper} applied to the i-th parameter. The implementation searches
    between \code{lower} and \code{upper} for the global optimum (minimum) of \code{fn}.}
  \item{control}{a list of control parameters; see \code{\link{DEoptim.control}}.}
  \item{...}{further arguments to be passed to \code{fn}.}
}
\details{
  \code{DEoptim} performs optimization (minimization) of \code{fn}. 

  The \code{control} argument is a list; see the help file for
  \code{\link{DEoptim.control}} for details.

  The \R implementation of Differential Evolution (DE), \pkg{DEoptim}, was first published on the Comprehensive \R Archive
  Network (CRAN) in 2005 by David Ardia. Early versions were written in 
  pure \R. Since version 2.0-0 (published to CRAN in 2009) the package has relied on an
  interface to a C implementation of DE, which is significantly 
  faster on most problems as compared to the implementation in 
  pure \R. The C interface is in many respects similar to the MS
  Visual C++ v5.0 implementation of the Differential Evolution algorithm
  distributed with the book \emph{Differential Evolution -- A Practical
  Approach to Global Optimization} by Price, K.V., Storn, R.M., Lampinen
  J.A, Springer-Verlag, 2006, and found on-line at
  \url{http://www.icsi.berkeley.edu/~storn/}. Since version 2.0-3 the C 
  implementation dynamically allocates the memory required to store the population, removing limitations on the 
  number of members in the population and length of the parameter vectors that may be optimized. 
  Since becoming publicly available, the package \pkg{DEoptim} has been used by several authors to solve optimization 
  problems arising in diverse domains; see Mullen et al. (2009) for a review.
  
  To perform a maximization (instead of minimization) of a given
  function, simply define a new function which is the opposite of the
  function to maximize and apply \code{DEoptim} to it.
  
  To integrate additional constraints (than box constraints) on the parameters \code{x} of
  \code{fn(x)}, for instance \code{x[1] + x[2]^2 < 2}, integrate the
  constraint within the function to optimize, for instance: 
  \preformatted{
    fn <- function(x)\{
      if (x[1] + x[2]^2 < 2)\{
        r <- Inf
      else\{
        ...
      \}
      return(r)
    \}
  }
  This simplistic strategy usually does not work all that well for gradient-based or Newton-type 
  methods. It is likely to be alright when the solution is in the interior of the feasible region, but when 
  the solution is on the boundary, optimization algorithm would have a difficult time converging. Furthermore, when
  the solution is on the boundary, this strategy would make the algorithm converge to an inferior solution in the interior. 
  However, for methods such as DE which are not gradient based, this strategy might not be that bad.

  Note that \code{DEoptim} stops if any \code{NA} or \code{NaN} value is
  obtained. You have to redefine your function to handle these values
  (for instance, set \code{NA} to \code{Inf} in your objective function).

  It is important to emphasize that the result of \code{DEoptim} is a random variable, 
  i.e., different results will obtain when the algorithm is run repeatedly with the same 
  settings. Hence, the user should set the random seed if they want to reproduce the results, e.g., by 
  setting \code{set.seed(1234)} before the call of \code{DEoptim}.

  \code{DEoptim} relies on repeated evaluation of the objective function
  in order to move the population toward a global minimum. Users
  interested in making \code{DEoptim} run as fast as possible should
  ensure that evaluation of the objective function is as efficient as
  possible. Using pure \R code, this may often be accomplished
  using vectorization. Writing parts of the objective function in a
  lower-level language like C or Fortran may also increase speed.

  Further details and examples of the \R package \pkg{DEoptim} can be found
  in Mullen et al. (2009) and Ardia et al. (2010).
  
  Please cite the package in publications. Use \code{citation("DEoptim")}.
}
\value{
  The output of the function \code{DEoptim} is a member of the \code{S3} class \code{DEoptim}. More precisely,
  this is a list (of length 2) containing the following elements:\cr

  \code{optim}, a list containing the following elements:
  \itemize{
    \item \code{bestmem}: the best set of parameters found.
    \item \code{bestval}: the value of \code{fn} corresponding to \code{bestmem}.
    \item \code{nfeval}: number of function evaluations.
    \item \code{iter}: number of procedure iterations.
  }

  \code{member}, a list containing the following elements:
  \itemize{
    \item \code{lower}: the lower boundary.
    \item \code{upper}: the upper boundary.
    \item \code{bestvalit}: the best value of \code{fn} at each iteration.
    \item \code{bestmemit}: the best member at each iteration.
    \item \code{pop}: the population generated at the last iteration.
    \item \code{storepop}: a list containing the intermediate populations.
  }

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
  \emph{Differential Evolution} (DE) is a search heuristic introduced by Storn and Price (1997). 
  Its remarkable performance as a global optimization algorithm on continuous numerical minimization 
  problems has been extensively explored; see Price et al. (2006). DE belongs to the class of genetic 
  algorithms which use biology-inspired operations of 
  crossover, mutation, and selection on a population in order to minimize an objective 
  function over the course of successive generations (see Mitchell, 1998). As with other evolutionary algorithms, 
  DE solves optimization problems by evolving a population of candidate solutions using alteration and selection
  operators. DE uses floating-point instead of bit-string encoding of population members, and 
  arithmetic operations instead of logical operations in mutation. DE is particularly well-suited to find the global optimum of a 
  real-valued function of real-valued parameters, and does not require that the function be
  either continuous or differentiable.

  Let \eqn{\mathit{NP}}{NP} denote the number of parameter vectors (members) \eqn{x \in R^d}{x in R^d} in the population. 
  In order to create the initial generation, \eqn{\mathit{NP}}{NP} guesses for the optimal value
  of the parameter vector are made, either using random values between lower and upper 
  bounds (defined by the user) or using values given by
  the user. Each generation involves creation of a new population from
  the current population members \eqn{\{ x_i \,|\, i = 1, \ldots, \mathit{NP}\}}{{x_i | i=1,...,NP}}, 
  where \eqn{i} indexes the vectors that make up the population.
  This is accomplished using \emph{differential mutation} of the
  population members. An initial mutant parameter vector \eqn{v_i} is
  created by choosing three members of the population, \eqn{x_{r_0}},
  \eqn{x_{r_1}} and \eqn{x_{r_2}}, at random. Then \eqn{v_i} is
  generated as

  \deqn{v_i \doteq x_{r_0} + \mathit{F} \cdot (x_{r_1} - x_{r_2})}{v_i := x_{r_0} + F * (x_{r_1} - x_{r_2})}  

  where \eqn{\mathit{F}}{F} is a positive scale factor, effective values for which are
  typically less than one. After the first mutation operation, mutation is
  continued until \eqn{d} mutations have been made, with a crossover probability 
  \eqn{\mathit{CR} \in [0,1]}{CR in [0,1]}.
  The crossover probability \eqn{\mathit{CR}}{CR} controls the fraction of the parameter
  values that are copied from the mutant. If an element of the trial parameter vector is found to violate the
  bounds after mutation and crossover, it is reset in such a way that the bounds are respected (with the
  specific protocol depending on the implementation).
  Then, the objective function values associated with the children are determined. If a trial
  vector has equal or lower objective function value than the previous
  vector it replaces the previous vector in the population;
  otherwise the previous vector remains. Variations of this scheme have also
  been proposed; see Price et al. (2006) and \code{\link{DEoptim.control}}.

  Intuitively, the effect of the scheme is that the shape of the distribution of the population in the
  search space is converging with respect to size and direction towards areas with high
  fitness. The closer the population gets to the global optimum, the more the distribution
  will shrink and therefore reinforce the generation of smaller difference vectors.  

  As a general advice regarding the choice of \eqn{\mathit{NP}}{NP}, \eqn{\mathit{F}}{F} and \eqn{\mathit{CR}}{CR}, 
  Storn et al. (2006) state the following: Set the number
  of parents \eqn{\mathit{NP}}{NP} to 10 times the number of parameters, select weighting factor 
  \eqn{\mathit{F} = 0.8}{F = 0.8} and crossover constant \eqn{\mathit{CR} = 0.9}{CR = 0.9}. Make sure that you initialize your parameter vectors
  by exploiting their full numerical range, i.e., if a parameter is allowed to exhibit
  values in the range [-100, 100] it is a good idea to pick the initial values from this
  range instead of unnecessarily restricting diversity. If you experience misconvergence in 
  the optimization process you usually have to increase the value for \eqn{\mathit{NP}}{NP}, but often you only have to adjust
  \eqn{\mathit{F}}{F} to be a little lower or higher than 0.8. If you increase
  \eqn{\mathit{NP}}{NP} and simultaneously lower \eqn{\mathit{F}}{F} a little, convergence is more
  likely to occur but generally takes longer, i.e., DE is getting
  more robust (there is always a convergence speed/robustness trade-off).

  DE is much more sensitive to the choice of \eqn{\mathit{F}}{F} than it is to
  the choice of \eqn{\mathit{CR}}{CR}. \eqn{\mathit{CR}}{CR} is more like a fine tuning element. High
  values of \eqn{\mathit{CR}}{CR} like \eqn{\mathit{CR} = 1}{CR = 1} give faster convergence if convergence
  occurs. Sometimes, however, you have to go down as much as \eqn{\mathit{CR} = 0}{CR = 0} to
  make DE robust enough for a particular problem. For more details on the DE strategy, we refer 
  the reader to Storn and Price (1997) and Price et al. (2006).
}
\references{
  Differential Evolution homepage:
  URL \url{http://www.icsi.berkeley.edu/~storn/code.html}
  
  Storn, R. and Price, K. (1997) 
  Differential Evolution -- A Simple and Efficient Heuristic for Global Optimization over Continuous Spaces,
  \emph{Journal of Global Optimization}, 11:4, 341--359.  

  Price, K.V., Storn, R.M., Lampinen J.A. (2006)
  \emph{Differential Evolution - A Practical Approach to Global Optimization}.
  Berlin Heidelberg: Springer-Verlag. ISBN 3540209506.

  Mitchell, M. (1998) 
  \emph{An Introduction to Genetic Algorithms}.
  The MIT Press. ISBN 0262631857.

  Mullen, K.M., Ardia, D., Gil, D.L, Windover, D., Cline, J. (2009)
  \pkg{DEoptim}: An \R Package for Global Optimization by Differential Evolution.
  URL \url{http://ssrn.com/abstract=1526466}

  Ardia, D., Boudt, K., Carl, P., Mullen, K.M., Peterson, B.G. (2010)
  Differential Evolution (\pkg{DEoptim}) for Non-Convex Portfolio Optimization.
  URL \url{http://ssrn.com/abstract=1584905}
}
\author{
  David Ardia, Katharine Mullen \email{katharine.mullen@nist.gov}, 
  Brian Peterson and Joshua Ulrich.
}
\seealso{
  \code{\link{DEoptim.control}} for control arguments,
  \code{\link{DEoptim-methods}} for methods on \code{DEoptim} objects,
  including some examples in plotting the results;
  \code{\link{optim}} or \code{\link{constrOptim}}
  for alternative optimization algorithms.
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

  ## DEoptim searches for minima of the objective function between
  ## lower and upper bounds on each parameter to be optimized. Therefore
  ## in the call to DEoptim we specify vectors that comprise the
  ## lower and upper bounds; these vectors are the same length as the
  ## parameter vector.
  lower <- c(-10,-10)
  upper <- -lower
 
  ## run DEoptim and set a seed first for replicability
  set.seed(1234)
  DEoptim(Rosenbrock, lower, upper)

  ## increase the population size
  DEoptim(Rosenbrock, lower, upper, DEoptim.control(NP = 100))

  ## change other settings and store the output
  outDEoptim <- DEoptim(Rosenbrock, lower, upper, DEoptim.control(NP = 80,
                        itermax = 400, F = 1.2, CR = 0.7))
  
  ## plot the output
  plot(outDEoptim)

  ## 'Wild' function, global minimum at about -15.81515
  Wild <- function(x)
    10 * sin(0.3 * x) * sin(1.3 * x^2) +
       0.00001 * x^4 + 0.2 * x + 80

  plot(Wild, -50, 50, n = 1000, main = "'Wild function'")

  outDEoptim <- DEoptim(Wild, lower = -50, upper = 50,
                        control = DEoptim.control(trace = FALSE))
  
  plot(outDEoptim)

  DEoptim(Wild, lower = -50, upper = 50,
          control = DEoptim.control(NP = 50))
}
\keyword{nonlinear}
\keyword{optimize}
