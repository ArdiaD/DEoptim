#include <R.h>
#include <Rinternals.h>

SEXP getListElement(SEXP list, char *str)
{
   SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
   int i;

   for (i = 0; i < length(list); i++)
     if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
       elmt = VECTOR_ELT(list, i);
       break;
     }
   return elmt;
}
