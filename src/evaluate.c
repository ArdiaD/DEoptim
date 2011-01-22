#include <R.h>
#include <Rdefines.h>

/*------objective function---------------------------------------*/

double evaluate(long *l_nfeval, SEXP par, SEXP fcall, SEXP env)
{
   SEXP sexp_fvec, fn;
   double f_result;

   PROTECT(fn = lang3(fcall, par, R_DotsSymbol));
      (*l_nfeval)++;  /* increment function evaluation count */

   PROTECT(sexp_fvec = eval(fn, env));
   f_result = NUMERIC_POINTER(sexp_fvec)[0];

   if(ISNAN(f_result))
     error("NaN value of objective function! \nPerhaps adjust the bounds.");

   UNPROTECT(2);
   return(f_result);
}
