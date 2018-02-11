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

typedef SEXP (*func_ptr)(SEXP, SEXP);
SEXP genrose(SEXP x_, SEXP e_)
{
  double *x = REAL(x_);
  int nr = nrows(x_);
  int nc = ncols(x_);

  SEXP s_a = findVar(install("a"), e_);
  SEXP s_b = findVar(install("b"), e_);

  if (s_a == R_UnboundValue) Rprintf("oops a\n");
  if (s_b == R_UnboundValue) Rprintf("oops b\n");

  double a = REAL(s_a)[0];
  double b = REAL(s_b)[0];

  SEXP s_out = PROTECT(allocVector(REALSXP, nr));
  double *out = REAL(s_out);

  double sum = 1.0;
  for (int i = 0; i < nr; i++) {
    sum = 1.0;
    int ij, ij1;
    for (int j = 1; j < nc; j++) {
      ij = i+nr*j;
      ij1 = i+nr*(j-1);
      sum += b*(pow(x[ij1]*x[ij1] - x[ij], 2)) + (x[ij] - a)*(x[ij] - a);
    }
    out[i] = sum; }

  UNPROTECT(1);
  return s_out;
}

SEXP popEvaluate(long *l_nfeval, SEXP parMat, SEXP fcall, SEXP env,
		 int incrementEval)
{
   SEXP sexp_fvec, fn;
   double *d_result;
   int P = 0;

   if (isNull(fcall))
     return parMat;

   if (TYPEOF(fcall) == EXTPTRSXP) {
     func_ptr fptr = NULL;
     fptr = (func_ptr)EXTPTR_PTR(fcall);
     PROTECT(sexp_fvec = fptr(parMat, env)); P++;
   } else {
     PROTECT(fn = lang3(fcall, parMat, R_DotsSymbol)); P++;
     PROTECT(sexp_fvec = eval(fn, env)); P++;
   }
   int nr = nrows(sexp_fvec);
   if(incrementEval)
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
