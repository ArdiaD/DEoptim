#include "de.h"
#include <R.h>
#include <Rdefines.h>

//------objective function---------------------------------------

t_pop evaluate(t_pop t_tmp, long *l_nfeval, t_pop *tpa_array, 
	            int i_NP, SEXP fcall, SEXP env)
{
   int i;
   SEXP par;
   SEXP sexp_fvec, fn;
   extern int gi_D; 
   float f_result;  ;

   PROTECT(par = NEW_NUMERIC(gi_D));
   for (i = 0; i < gi_D; i++) {
      NUMERIC_POINTER(par)[i] = (double)t_tmp.fa_vector[i];
   }
   PROTECT(fn = lang2(fcall, par));
      (*l_nfeval)++;  //increment function evaluation count
  
   SETCADR(fn, par);
 
   PROTECT(sexp_fvec = eval(fn, env));
   f_result = (float)(NUMERIC_POINTER(sexp_fvec)[0]);
 
   UNPROTECT(3);	
 
   t_tmp.fa_cost[0] = f_result;
 
   return(t_tmp);
}

int left_vector_wins(t_pop t_trial, t_pop t_target)
{
   //---trial wins against target even when cost is equal.-----
   if (t_trial.fa_cost[0] <= t_target.fa_cost[0]) 
      return(TRUE);
   else 
      return(FALSE);
}
