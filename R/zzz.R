".onLoad" <- function (lib, pack)
{
  library.dynam(pack, pack, lib)
  packageStartupMessage("\nDEoptim package",
   "\nDifferential Evolution algorithm in R",
   "\nAuthors: David Ardia and Katharine Mullen\n")
}
