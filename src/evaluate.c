#include <R.h>
#include <Rdefines.h>

/*------objective function---------------------------------------*/

double evaluate(long *l_nfeval, double *param, SEXP par, SEXP fcall, SEXP env)
{
   int i, P;
   SEXP sexp_fvec, fn;
   double f_result;  

   for (i = 0; i < nrows(par); i++) {
      NUMERIC_POINTER(par)[i] = param[i];
   }
   PROTECT(fn = lang2(fcall, par)); P++;
      (*l_nfeval)++;  /* increment function evaluation count */

   PROTECT(sexp_fvec = eval(fn, env)); P++;
   f_result = NUMERIC_POINTER(sexp_fvec)[0];

   if(ISNAN(f_result))
     error("NaN value of objective function! \nPerhaps adjust the bounds.");

   UNPROTECT(P);

   return(f_result); 
}
