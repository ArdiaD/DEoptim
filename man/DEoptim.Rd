\name{DEoptim}
\alias{DEoptim}
\title{Differential Evolution Optimization}
\concept{minimization}
\description{
  Performs evolutionary global optimization via the differential evolution algorithm.
}
\usage{
DEoptim(fn, lower, upper, control = DEoptim.control(), ...)
}
\arguments{
  \item{fn}{A function to be minimized, with first argument the vector
    of parameters over which minimization is to take place. It should
    return a scalar result. \code{NA} and \code{NaN} values are not
    allowed.}
  \item{lower, upper}{Bounds on the variables.}
  \item{control}{A list of control parameters;
    see \code{\link{DEoptim.control}}.}
  \item{...}{Further arguments to be passed to \code{fn}.}
}
\details{
  \code{DEoptim} performs minimization of \code{fn}. 

  The \code{control} argument is a list; see the help file for
  \code{\link{DEoptim.control}} for details.

}
\value{
  A list of lists of the class \code{DEoptim}.\cr

  list \code{optim} contains the followings:\cr
  \code{bestmem}: the best set of parameters found.\cr
  \code{bestval}: the value of \code{fn} corresponding to \code{bestmem}.\cr
  \code{nfeval}: number of function evaluations.\cr
  \code{iter}: number of procedure iterations.\cr

  list \code{member} contains the followings:\cr
  \code{lower}: the lower boundary.\cr
  \code{upper}: the upper boundary.\cr
  \code{bestvalit}: the best value of \code{fn} at each iteration.\cr
  \code{bestmemit}: the best member at each iteration.\cr
  \code{pop}: the population generated at the last iteration.\cr
  \code{storepop}: a list containing the intermediate populations.\cr
}
\note{
  \code{DEoptim} calls C code that is in many respects similar to the MS
  Visual C++ v5.0 implementation of the differential evolution algorithm
  distributed with the book \emph{Differential Evolution -- A Practical
  Approach to Global Optimization} by Price, K.V., Storn, R.M., Lampinen
  J.A, Springer-Verlag, 2005, and found on-line at
  \url{http://www.icsi.berkeley.edu/~storn/DeWin.zip}.

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
  
  To perform a maximization (instead of minimization) of a given
  function, simply define a new function which is the opposite of the
  function to maximize and apply \code{DEoptim} to it.
  
  To integrate additional constraints on the parameters \code{x} of
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

  Note that \code{DEoptim} stops if any \code{NA} or \code{NaN} value is
  obtained. You have to redefine your function to handle these values
  (for instance, set \code{NA} to \code{Inf} in your objective function).
  
  You can reproduce your results by setting the seed of the random number 
  generator, e.g. \code{set.seed(1234)}.
  
  Please cite the package in publications. Use \code{citation("DEoptim")}.
}
\references{
  Differential evolution homepage:
  \url{http://www.icsi.berkeley.edu/~storn/code.html}
  
  Storn, R. and Price, K. (1997). Differential Evolution -- A Simple and
  Efficient Heuristic for Global Optimization over Continuous Spaces,
  \emph{Journal of Global Optimization}, 11:4, 341-359.  

  Price, K.V., Storn, R.M., Lampinen J.A. (2005).
  \emph{Differential Evolution - A Practical Approach to Global Optimization}.
  Springer-Verlag. ISBN 3540209506.

  Mullen, K.M., Ardia, D., Gil, D.L, Windover, D., Cline, J. (2009).
  DEoptim: An R Package for Global Optimization by Differential Evolution.
  \url{http://ssrn.com/abstract=1526466}

  Ardia, D., Boudt, K., Carl, P., Mullen, K.M., Peterson, B.G. (2010).
  Differential Evolution (DEoptim) for Non-Convex Portfolio Optimization.
  \url{http://ssrn.com/abstract=1584905}
}
\author{
  David Ardia \email{david.ardia@unifr.ch} and 
  Katharine Mullen \email{katharine.mullen@nist.gov}.
}
\note{%
  Further details and examples of the \R package \code{DEoptim} can be found
  in Mullen et al. (2009) and Ardia et al. (2010).
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
  Rosenbrock <- function(x){
    x1 <- x[1]
    x2 <- x[2]
    100 * (x2 - x1 * x1)^2 + (1 - x1)^2
  }

  lower <- c(-10,-10)
  upper <- -lower
 
  set.seed(1234)

  DEoptim(Rosenbrock, lower, upper)

  DEoptim(Rosenbrock, lower, upper, DEoptim.control(NP = 100))

  outDEoptim <- DEoptim(Rosenbrock, lower, upper, DEoptim.control(NP = 80,
                        itermax = 400, F = 1.2, CR = 0.7))
  
  plot(outDEoptim)

  ## 'Wild' function, global minimum at about -15.81515
  Wild <- function(x)
    10 * sin(0.3 * x) * sin(1.3 * x^2) +
       0.00001 * x^4 + 0.2 * x + 80

  
  plot(Wild, -50, 50, n = 1000,
    main = "'Wild function'")

  outDEoptim <- DEoptim(Wild, lower = -50, upper = 50,
                        control = DEoptim.control(trace = FALSE))
  
  plot(outDEoptim)

  DEoptim(Wild, lower = -50, upper = 50,
          control = DEoptim.control(NP = 50))
}
\keyword{nonlinear}
\keyword{optimize}