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
   int P = 0;

   PROTECT(fn = lang3(fcall, parMat, R_DotsSymbol)); P++;
 
   PROTECT(sexp_fvec = eval(fn, env)); P++;
   int nr = nrows(sexp_fvec);
   (*l_nfeval) += nr;  
   if(nr != nrows(parMat))
     error("objective function result has different length than parameter matrix");
   switch(TYPEOF(sexp_fvec)) {
     case INTSXP:
       PROTECT(sexp_fvec = coerceVector(sexp_fvec, REALSXP)); P++;
       break;
     case REALSXP:
       break;
     default:
       error("unsupported objective function return value");
       break;
   }
   d_result = REAL(sexp_fvec);
   for(int i=0; i < nr; i++) {
     if(ISNAN(d_result[i]))
       error("NaN value of objective function! \nPerhaps adjust the bounds.");
   }
   UNPROTECT(P);
   return(sexp_fvec);
}
