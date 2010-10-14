#include <R.h>
#include <Rdefines.h>

/*------objective function---------------------------------------*/

double evaluate(long *l_nfeval, double *param, SEXP par, SEXP fcall, SEXP env)
{
   int i;
   SEXP sexp_fvec, fn;
   double f_result;  

   for (i = 0; i < nrows(par); i++) {
      NUMERIC_POINTER(par)[i] = param[i];
   }
   PROTECT(fn = lang2(fcall, par)); 
      (*l_nfeval)++;  /* increment function evaluation count */

   PROTECT(sexp_fvec = eval(fn, env)); 
   f_result = NUMERIC_POINTER(sexp_fvec)[0];
   
   UNPROTECT(2);
   if(ISNAN(f_result))
     error("NaN value of objective function! \nPerhaps adjust the bounds.");
   
   return(f_result); 
}
