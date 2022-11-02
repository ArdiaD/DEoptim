library(DEoptim)

# Tests for partial argument matching (#2)
# 'p' matched 'params' (prior to ae82a815751f6c1eeaef4e51430df8fcef556af9)
# 'M' matched MARGIN
# 'F' matched FUN
RosenbrockPartialMatch <-
function(x, p, M, F)
{
  x1 <- x[1]
  x2 <- x[2]
  100 * (x2 - x1 * x1)^2 + (1 - x1)^2 
}
lower <- c(-10, -10)
upper <- -lower

set.seed(21)
deCtrl <- DEoptim.control(trace = FALSE), parallelType = "parallel")
deopt <- DEoptim(RosenbrockPartialMatch, lower, upper, deCtrl,
                 p = NULL, M = NULL, F = NULL)