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

SEXP popEvaluate(long *l_nfeval, SEXP parMat, SEXP fcall, SEXP env)
{
   SEXP sexp_fvec, fn;
   double *d_result;
   int *i_result;

   PROTECT(fn = lang3(fcall, parMat, R_DotsSymbol));
      (*l_nfeval)++;  /* increment function evaluation count */

   PROTECT(sexp_fvec = eval(fn, env));
   int nr = nrows(sexp_fvec);
   switch(TYPEOF(sexp_fvec)) {
     case INTSXP:
       i_result = INTEGER(sexp_fvec);
       for(int i=0; i < nr; i++) {
         if(ISNAN(i_result[i]))
           error("NaN value of objective function! \nPerhaps adjust the bounds.");
       }
       break;
     case REALSXP:
       d_result = REAL(sexp_fvec);
       for(int i=0; i < nr; i++) {
         if(ISNAN(d_result[i]))
           error("NaN value of objective function! \nPerhaps adjust the bounds.");
       }
       break;
     default:
       error("unsupported objective function return value");
       break;
   }
   UNPROTECT(2);
   return(sexp_fvec);
}
