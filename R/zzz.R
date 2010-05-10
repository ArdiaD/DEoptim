".onLoad" <- function (lib, pack)
{
  library.dynam(pack, pack, lib)
  cat("\nDEoptim package")
  cat("\nDifferential Evolution algorithm in R")
  cat("\nAuthors: David Ardia and Katharine Mullen\n")
}
