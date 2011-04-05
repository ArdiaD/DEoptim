".onLoad" <- function (lib, pack)
{
  library.dynam(pack, pack, lib)
  packageStartupMessage("\nDEoptim package",
   "\nDifferential Evolution algorithm in R",
   "\nAuthors: D. Ardia, K. Mullen, B. Peterson and J. Ulrich\n")
}
