#include <R.h>
#include <Rdefines.h>

/*------objective function---------------------------------------*/

//double evaluate(long *l_nfeval, double *param, int i_D, SEXP fcall, SEXP env)
double evaluate(double *param, int i_D, SEXP fcall, SEXP env)
{
   int i;
   SEXP par;
   SEXP sexp_fvec, fn;
   double f_result;  

   PROTECT(par = NEW_NUMERIC(i_D));
   for (i = 0; i < i_D; i++) {
      NUMERIC_POINTER(par)[i] = param[i];
   }
   PROTECT(fn = lang2(fcall, par));
      //(*l_nfeval)++;  /* increment function evaluation count */
  
   SETCADR(fn, par);
 
   PROTECT(sexp_fvec = eval(fn, env));
   f_result = NUMERIC_POINTER(sexp_fvec)[0];
 
   UNPROTECT(3);	
  
   if(ISNAN(f_result))
     error("NaN value of objective function! \nPerhaps adjust the bounds.");
   
   return(f_result); 
}
