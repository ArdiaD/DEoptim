 /***************************************************************

Implementation of DE based loosely on DE-Engine v4.0, Rainer Storn, 2004
by Katharine Mullen, 2009
modified by Joshua Ulrich, 2010

Storn's MS Visual C++ v5.0 was downloaded from
http://www.icsi.berkeley.edu/~storn/DeWin.zip

Note that the struct t_pop was not used since we want to store in it a
parameter vector of arbitrary length --

using a construction like
typedef struct t_pop
{
  double fa_cost[MAXCOST];
  double fa_vector[];
} t_pop;
I have not been able to use Calloc or R_alloc to set aside memory for
type t_pop dynamically; I did not look into the bowels but believe this
is likely a limitation of these functions.

***************************************************************/

/*------Include section----------------------------------------------*/

#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>

SEXP getListElement(SEXP list, char *str);
SEXP DEoptimC(SEXP lower, SEXP upper, SEXP fn, SEXP control, SEXP rho, SEXP fnMap);
void devol(double VTR, double d_weight, double fcross, int i_bs_flag,
           double *d_lower, double *d_upper, SEXP fcall, SEXP rho, int i_trace,
           int i_strategy, int i_D, int i_NP, int i_itermax,
           double *initialpopv, int i_storepopfreq, int i_storepopfrom,
           int i_specinitialpop, 
           double *gt_bestP, double *gt_bestC,
           double *gd_pop, double *gd_storepop, double *gd_bestmemit, double *gd_bestvalit,
           int *gi_iter, double d_pPct, double d_c, long *l_nfeval,
           double d_reltol, int i_steptol, SEXP fnMap);
void permute(int ia_urn2[], int i_urn2_depth, int i_NP, int i_avoid, int ia_urn1[]);
double evaluate(long *l_nfeval, SEXP par, SEXP fcall, SEXP env);
SEXP popEvaluate(long *l_nfeval, SEXP parMat, SEXP fcall, SEXP env, int incrementEval);


/*------General functions-----------------------------------------*/

SEXP DEoptimC(SEXP lower, SEXP upper, SEXP fn, SEXP control, SEXP rho, SEXP fnMap)
{
  int i, j, P=0;

  if (!isFunction(fn))
    error("fn is not a function!");
  if (!isEnvironment(rho))
    error("rho is not an environment!");

  /*-----Initialization of annealing parameters-------------------------*/
  /* value to reach */
  double VTR = NUMERIC_VALUE(getListElement(control, "VTR"));
  /* chooses DE-strategy */
  int i_strategy = INTEGER_VALUE(getListElement(control, "strategy"));
  /* Maximum number of generations */
  int i_itermax = INTEGER_VALUE(getListElement(control, "itermax"));
  /* Dimension of parameter vector */
  int i_D = INTEGER_VALUE(getListElement(control, "npar"));
  /* Number of population members */
  int i_NP = INTEGER_VALUE(getListElement(control, "NP"));
  /* When to start storing populations */
  int i_storepopfrom = INTEGER_VALUE(getListElement(control, "storepopfrom"))-1;
  /* How often to store populations */
  int i_storepopfreq = INTEGER_VALUE(getListElement(control, "storepopfreq"));
  /* User-defined inital population */
  int i_specinitialpop = INTEGER_VALUE(getListElement(control, "specinitialpop"));
  double *initialpopv = NUMERIC_POINTER(getListElement(control, "initialpop"));
  /* stepsize */
  double d_weight = NUMERIC_VALUE(getListElement(control, "F"));
  /* crossover probability */
  double d_cross = NUMERIC_VALUE(getListElement(control, "CR"));
  /* Best of parent and child */
  int i_bs_flag = NUMERIC_VALUE(getListElement(control, "bs"));
  /* Print progress? */
  int i_trace = NUMERIC_VALUE(getListElement(control, "trace"));
  /* p to define the top 100p% best solutions */
  double d_pPct = NUMERIC_VALUE(getListElement(control, "p"));
  /* crossover adaptation (a positive constant between 0 and 1) */
  double d_c = NUMERIC_VALUE(getListElement(control, "c"));
  /* relative tolerance */
  double d_reltol = NUMERIC_VALUE(getListElement(control, "reltol"));
  /* relative tolerance steps */
  int i_steptol = NUMERIC_VALUE(getListElement(control, "steptol"));

  int i_nstorepop = ceil((i_itermax - i_storepopfrom) / i_storepopfreq);
  /* Use S_alloc, since it initializes with zeros FIXME: these should be SEXP */
  double *gd_storepop = (double *)S_alloc(i_NP,sizeof(double) * i_D * i_nstorepop);

  /* External pointers to return to R */
  SEXP sexp_bestmem, sexp_bestval, sexp_nfeval, sexp_iter,
    out, sexp_pop, sexp_storepop, sexp_bestmemit, sexp_bestvalit;

  PROTECT(sexp_bestmem = NEW_NUMERIC(i_D)); P++;
  PROTECT(sexp_pop = allocMatrix(REALSXP, i_D, i_NP)); P++;
  PROTECT(sexp_bestmemit = allocMatrix(REALSXP, i_itermax, i_D)); P++;
  PROTECT(sexp_bestvalit = allocVector(REALSXP, i_itermax)); P++;
  double *gt_bestP     = REAL(sexp_bestmem);
  double *gd_pop       = REAL(sexp_pop);
  double *gd_bestmemit = REAL(sexp_bestmemit);
  double *gd_bestvalit = REAL(sexp_bestvalit);

  /* ensure lower and upper are double */
  if(TYPEOF(lower) != REALSXP) {PROTECT(lower = coerceVector(lower, REALSXP)); P++;}
  if(TYPEOF(upper) != REALSXP) {PROTECT(upper = coerceVector(upper, REALSXP)); P++;}
  double *d_lower      = REAL(lower);
  double *d_upper      = REAL(upper);

  double gt_bestC;
  int gi_iter = 0;
  long l_nfeval = 0;

  /*---optimization--------------------------------------*/
  devol(VTR, d_weight, d_cross, i_bs_flag, d_lower, d_upper, fn, rho, i_trace,
        i_strategy, i_D, i_NP, i_itermax,
        initialpopv, i_storepopfrom, i_storepopfreq,
        i_specinitialpop,
        gt_bestP, &gt_bestC,
        gd_pop, gd_storepop, gd_bestmemit, gd_bestvalit,
        &gi_iter, d_pPct, d_c, &l_nfeval,
        d_reltol, i_steptol, fnMap);
  /*---end optimization----------------------------------*/

  j =  i_nstorepop * i_NP * i_D;
  PROTECT(sexp_storepop = NEW_NUMERIC(j)); P++;
  for (i = 0; i < j; i++)
    NUMERIC_POINTER(sexp_storepop)[i] = gd_storepop[i];

  PROTECT(sexp_nfeval = ScalarInteger(l_nfeval)); P++;
  PROTECT(sexp_iter = ScalarInteger(gi_iter)); P++;
  PROTECT(sexp_bestval = ScalarReal(gt_bestC)); P++;

  const char *out_names[] = {"bestmem", "bestval", "nfeval",
      "iter", "bestmemit", "bestvalit", "pop", "storepop", ""};
  PROTECT(out = mkNamed(VECSXP, out_names)); P++;
  SET_VECTOR_ELT(out, 0, sexp_bestmem);
  SET_VECTOR_ELT(out, 1, sexp_bestval);
  SET_VECTOR_ELT(out, 2, sexp_nfeval);
  SET_VECTOR_ELT(out, 3, sexp_iter);
  SET_VECTOR_ELT(out, 4, sexp_bestmemit);
  SET_VECTOR_ELT(out, 5, sexp_bestvalit);
  SET_VECTOR_ELT(out, 6, sexp_pop);
  SET_VECTOR_ELT(out, 7, sexp_storepop);

  UNPROTECT(P);
  return out;
}

void devol(double VTR, double d_weight, double d_cross, int i_bs_flag,
           double *d_lower, double *d_upper, SEXP fcall, SEXP rho, int trace,
           int i_strategy, int i_D, int i_NP, int i_itermax,
           double *initialpopv, int i_storepopfrom, int i_storepopfreq,
           int i_specinitialpop, 
           double *gt_bestP, double *gt_bestC,
           double *gd_pop, double *gd_storepop, double *gd_bestmemit, double *gd_bestvalit,
           int *gi_iter, double d_pPct, double d_c, long *l_nfeval,
           double d_reltol, int i_steptol, SEXP fnMap)
{

#define URN_DEPTH  5   /* 4 + one index to avoid */

  int P=0;
  /* initialize parameter vector to pass to evaluate function */
  SEXP par;
  PROTECT(par = NEW_NUMERIC(i_D)); P++;
  
  /* Data structures for parameter vectors */
  SEXP sexp_gta_popP, sexp_gta_oldP, sexp_gta_newP, sexp_map_pop;
  PROTECT(sexp_gta_popP = allocMatrix(REALSXP, i_NP, i_D)); P++; /* FIXME THIS HAD 2x the rows!!! */
  PROTECT(sexp_gta_oldP = allocMatrix(REALSXP, i_NP, i_D)); P++;
  PROTECT(sexp_gta_newP = allocMatrix(REALSXP, i_NP, i_D)); P++;
  double *ngta_popP = REAL(sexp_gta_popP); /* FIXME THIS HAD 2x the rows!!! */
  double *ngta_oldP = REAL(sexp_gta_oldP);
  double *ngta_newP = REAL(sexp_gta_newP);

  /* Data structures for objective function values associated with
   * parameter vectors */
  SEXP sexp_gta_popC, sexp_gta_oldC, sexp_gta_newC;
  PROTECT(sexp_gta_popC = allocVector(REALSXP, i_NP)); P++;
  PROTECT(sexp_gta_oldC = allocVector(REALSXP, i_NP)); P++;
  PROTECT(sexp_gta_newC = allocVector(REALSXP, i_NP)); P++;
  double *ngta_popC = REAL(sexp_gta_popC);
  double *ngta_oldC = REAL(sexp_gta_oldC);
  double *ngta_newC = REAL(sexp_gta_newC);

  double *t_bestitP = (double *)R_alloc(1,sizeof(double) * i_D);
  
  SEXP sexp_t_tmpP, sexp_t_tmpC;
  PROTECT(sexp_t_tmpP = allocMatrix(REALSXP, i_NP, i_D)); P++;
  PROTECT(sexp_t_tmpC = allocVector(REALSXP, i_NP)); P++;
  double *nt_tmpP = REAL(sexp_t_tmpP);
  double *nt_tmpC = REAL(sexp_t_tmpC);

  int i, j, k;  /* counting variables */
  int i_r1, i_r2, i_r3;  /* placeholders for random indexes */

  int ia_urn2[URN_DEPTH];
  int ia_urnTemp[i_NP];
  
  int popcnt, bestacnt; /* lazy cnters */

  double d_jitter, d_dither;

  double t_bestC;

  double **initialpop = (double **)R_alloc(i_NP,sizeof(double *));
  for (int i = 0; i < i_NP; i++)
    initialpop[i] = (double *)R_alloc(i_D,sizeof(double));

  /* vars for DE/current-to-p-best/1 */
  int i_pbest;
  int p_NP = round(d_pPct * i_NP);  /* choose at least two best solutions */
      p_NP = p_NP < 2 ? 2 : p_NP;
  int sortIndex[i_NP];              /* sorted values of gta_oldC */
  for(i = 0; i < i_NP; i++) sortIndex[i] = i;
  //double goodCR = 0, goodF = 0, goodF2 = 0, meanCR = 0.5, meanF = 0.5;
  double goodCR = 0, goodF = 0, goodF2 = 0, meanCR = d_cross, meanF = d_weight;
  int i_goodNP = 0;

  /* vars for when i_bs_flag == 1 */
//  int i_len, done, step, bound;
//  double tempC;

  GetRNGstate();

  /* if initial population provided, initialize with values */
  if (i_specinitialpop > 0) {
    k = 0;

    for (j = 0; j < i_D; j++) {
      for (i = 0; i < i_NP; i++) {
        initialpop[i][j] = initialpopv[k];
        k += 1;
      }
    }
  }

  /*------Initialization-----------------------------*/
  for (j = 0; j < i_D; j++) {
    for (i = 0; i < i_NP; i++) {
      if (i_specinitialpop <= 0) { /* random initial member */
        ngta_popP[i+i_NP*j] = d_lower[j] +
        unif_rand() * (d_upper[j] - d_lower[j]);

      }
      else /* or user-specified initial member */
        ngta_popP[i+i_NP*j] = initialpop[i][j];
    }
  }
  PROTECT(sexp_map_pop  = popEvaluate(l_nfeval, sexp_gta_popP, fnMap, rho, 0));
  memmove(REAL(sexp_gta_popP), REAL(sexp_map_pop), i_NP * i_D * sizeof(double)); // valgrind reports memory overlap here
  UNPROTECT(1);  // sexp_map_pop
  PROTECT(sexp_gta_popC = popEvaluate(l_nfeval, sexp_gta_popP,  fcall, rho, 1));
  ngta_popC = REAL(sexp_gta_popC);
  for (i = 0; i < i_NP; i++) {
    if (i == 0 || ngta_popC[i] <= t_bestC) {
      t_bestC = ngta_popC[i];
      for (j = 0; j < i_D; j++)
        gt_bestP[j]=ngta_popP[i+i_NP*j];
    }
  }

  /*---assign pointers to current ("old") population---*/
  memmove(REAL(sexp_gta_oldP), REAL(sexp_gta_popP), i_NP * i_D * sizeof(double));
  memmove(REAL(sexp_gta_oldC), REAL(sexp_gta_popC), i_NP * sizeof(double));
  UNPROTECT(1);  // sexp_gta_popC

  /*------Iteration loop--------------------------------------------*/
  int i_iter = 0;
  popcnt = 0;
  bestacnt = 0;
  int i_iter_tol = 0;

  while ((i_iter < i_itermax) && (t_bestC > VTR) && (i_iter_tol <= i_steptol))
  {
    /* store intermediate populations */
    if (i_iter % i_storepopfreq == 0 && i_iter >= i_storepopfrom) {
      for (i = 0; i < i_NP; i++) {
        for (j = 0; j < i_D; j++) {
          gd_storepop[popcnt] = ngta_oldP[i+i_NP*j];
          popcnt++;
        }
      }
    } /* end store pop */

    /* store the best member */
    for(j = 0; j < i_D; j++) {
      gd_bestmemit[bestacnt] = gt_bestP[j];
      bestacnt++;
    }
    /* store the best value */
    gd_bestvalit[i_iter] = t_bestC;

    for (j = 0; j < i_D; j++)
      t_bestitP[j] = gt_bestP[j];

    i_iter++;

    /*----compute dithering factor -----------------*/
    if (i_strategy == 5)
      d_dither = d_weight + unif_rand() * (1.0 - d_weight);

    /*---DE/current-to-p-best/1 ----------------------------------------------*/
    if (i_strategy == 6) {
      /* create a copy of gta_oldC to avoid changing it */
      double temp_oldC[i_NP];
      for(j = 0; j < i_NP; j++) temp_oldC[j] = ngta_oldC[j];

      /* sort temp_oldC to use sortIndex later */
      rsort_with_index( (double*)temp_oldC, (int*)sortIndex, i_NP );
    }

    /*----start of loop through ensemble------------------------*/
    for (i = 0; i < i_NP; i++) {

      /*t_tmpP is the vector to mutate and eventually select*/
      for (j = 0; j < i_D; j++)
        nt_tmpP[i+i_NP*j] = ngta_oldP[i+i_NP*j];
      nt_tmpC[i] = ngta_oldC[i];

      permute(ia_urn2, URN_DEPTH, i_NP, i, ia_urnTemp); /* Pick 4 random and distinct */

      i_r1 = ia_urn2[1];  /* population members */
      i_r2 = ia_urn2[2];
      i_r3 = ia_urn2[3];

      if (d_c > 0) {
        d_cross = rnorm(meanCR, 0.1);
        d_cross = d_cross > 1.0 ? 1 : d_cross;
        d_cross = d_cross < 0.0 ? 0 : d_cross;
        do {
          d_weight = rcauchy(meanF, 0.1);
          d_weight = d_weight > 1 ? 1.0 : d_weight;
        }while(d_weight <= 0.0);
      }

      /*===Choice of strategy===============================================*/
      j = (int)(unif_rand() * i_D); /* random parameter */
      k = 0;
      do {
        switch (i_strategy) {
          case 1: { /*---classical strategy DE/rand/1/bin-------------------*/
            nt_tmpP[i+i_NP*j] = ngta_oldP[i_r1+i_NP*j] +
              d_weight * (ngta_oldP[i_r2+i_NP*j] - ngta_oldP[i_r3+i_NP*j]);
            break;
          }
          case 2: { /*---DE/local-to-best/1/bin-----------------------------*/
            nt_tmpP[i+i_NP*j] = nt_tmpP[i+i_NP*j] +
              d_weight * (t_bestitP[j] - nt_tmpP[i+i_NP*j]) +
              d_weight * (ngta_oldP[i_r2+i_NP*j] - ngta_oldP[i_r3+i_NP*j]);
            break;
          }
          case 3: { /*---DE/best/1/bin with jitter--------------------------*/
            d_jitter = 0.0001 * unif_rand() + d_weight;
            nt_tmpP[i+i_NP*j] = t_bestitP[j] +
              d_jitter * (ngta_oldP[i_r1+i_NP*j] - ngta_oldP[i_r2+i_NP*j]);
            break;
          }
          case 4: { /*---DE/rand/1/bin with per-vector-dither---------------*/
            nt_tmpP[i+i_NP*j] = ngta_oldP[i_r1+i_NP*j] +
              (d_weight + unif_rand()*(1.0 - d_weight))*
              (ngta_oldP[i_r2+i_NP*j]-ngta_oldP[i_r3+i_NP*j]);
            break;
          }
          case 5: { /*---DE/rand/1/bin with per-generation-dither-----------*/
            nt_tmpP[i+i_NP*j] = ngta_oldP[i_r1+i_NP*j] +
              d_dither * (ngta_oldP[i_r2+i_NP*j] - ngta_oldP[i_r3+i_NP*j]);
            break;
          }
          case 6: { /*---DE/current-to-p-best/1 (JADE)----------------------*/
            /* select from [0, 1, 2, ..., (pNP-1)] */
            i_pbest = sortIndex[(int)(unif_rand() * p_NP)];
            nt_tmpP[i+i_NP*j] = ngta_oldP[i+i_NP*j] +
              d_weight * (ngta_oldP[i_pbest+i_NP*j] - ngta_oldP[i+i_NP*j]) +
              d_weight * (ngta_oldP[i_r1+i_NP*j]    - ngta_oldP[i_r2+i_NP*j]);
            break;
          }
          default: { /*---variation to DE/rand/1/bin: either-or-algorithm---*/
            if (unif_rand() < 0.5) { /* differential mutation, Pmu = 0.5 */
              nt_tmpP[i+i_NP*j] = ngta_oldP[i_r1+i_NP*j] + d_weight *
                (ngta_oldP[i_r2+i_NP*j] - ngta_oldP[i_r3+i_NP*j]);
            } else {
            /* recombination with K = 0.5*(F+1) -. F-K-Rule */
            nt_tmpP[i+i_NP*j] = ngta_oldP[i_r1+i_NP*j] +
              0.5 * (d_weight + 1.0) * (ngta_oldP[i_r2+i_NP*j]
              + ngta_oldP[i_r3+i_NP*j] - 2 * ngta_oldP[i_r1+i_NP*j]);
            }
          }
        } /* end switch */
        j = (j + 1) % i_D;
        k++;
      }while((unif_rand() < d_cross) && (k < i_D));
      /*===End choice of strategy===========================================*/

      /*----boundary constraints, bounce-back method was not enforcing bounds correctly*/
      for (j = 0; j < i_D; j++) {
        if (nt_tmpP[i+i_NP*j] < d_lower[j]) {
          nt_tmpP[i+i_NP*j] = d_lower[j] + unif_rand() * (d_upper[j] - d_lower[j]);
        }
        if (nt_tmpP[i+i_NP*j] > d_upper[j]) {
          nt_tmpP[i+i_NP*j] = d_upper[j] - unif_rand() * (d_upper[j] - d_lower[j]);
        }
      }

    } /* NEW End mutation loop through ensemble */

    /*------Trial mutation now in t_tmpP-----------------*/
    /* evaluate mutated population */
    if(i_iter > 1) UNPROTECT(1);  // previous iteration's sexp_t_tmpC
    PROTECT(sexp_map_pop = popEvaluate(l_nfeval, sexp_t_tmpP,  fnMap, rho, 0));
    memmove(REAL(sexp_t_tmpP), REAL(sexp_map_pop), i_NP * i_D * sizeof(double)); // valgrind reports memory overlap here
    UNPROTECT(1);  // sexp_map_pop
    PROTECT(sexp_t_tmpC  = popEvaluate(l_nfeval, sexp_t_tmpP, fcall, rho, 1));
    nt_tmpC = REAL(sexp_t_tmpC);

    /* compare old pop with mutated pop */
    for (i = 0; i < i_NP; i++) {

      /* note that i_bs_flag means that we will choose the
       *best NP vectors from the old and new population later*/
      if (nt_tmpC[i] <= ngta_oldC[i] || i_bs_flag) {
        /* replace target with mutant */
        for (j = 0; j < i_D; j++)
          ngta_newP[i+i_NP*j]=nt_tmpP[i+i_NP*j];
        ngta_newC[i]=nt_tmpC[i];
        if (nt_tmpC[i] <= t_bestC) {
          for (j = 0; j < i_D; j++)
            gt_bestP[j]=nt_tmpP[i+i_NP*j];
          t_bestC=nt_tmpC[i];
        }
        if (d_c > 0) { /* calculate new goodCR and goodF */
          goodCR += d_cross / ++i_goodNP;
          goodF += d_weight;
          goodF2 += pow(d_weight,2.0);
        }
      }
      else {
        for (j = 0; j < i_D; j++)
          ngta_newP[i+i_NP*j]=ngta_oldP[i+i_NP*j];
        ngta_newC[i]=ngta_oldC[i];

      }
    } /* End mutation loop through ensemble */
   
    if (d_c > 0) { /* calculate new meanCR and meanF */
      meanCR = (1-d_c)*meanCR + d_c*goodCR;
      meanF = (1-d_c)*meanF + d_c*goodF2/goodF;
    }

    if(i_bs_flag) {  /* FIXME */
      error("bs = TRUE not currently supported");
//      /* examine old and new pop. and take the best NP members
//       * into next generation */
//      for (i = 0; i < i_NP; i++) {
//        for (j = 0; j < i_D; j++)
//          gta_popP[i][j] = gta_oldP[i][j];
//        gta_popC[i] = gta_oldC[i];
//      }
//      for (i = 0; i < i_NP; i++) {
//        for (j = 0; j < i_D; j++)
//          gta_popP[i_NP+i][j] = gta_newP[i][j];
//        gta_popC[i_NP+i] = gta_newC[i];
//      }
//      i_len = 2 * i_NP;
//      step = i_len;  /* array length */
//      while (step > 1) {
//        step /= 2;   /* halve the step size */
//        do {
//          done = 1;
//          bound  = i_len - step;
//          for (j = 0; j < bound; j++) {
//              i = j + step + 1;
//              if (gta_popC[j] > gta_popC[i-1]) {
//                  for (k = 0; k < i_D; k++)
//                    tempP[k] = gta_popP[i-1][k];
//                  tempC = gta_popC[i-1];
//                  for (k = 0; k < i_D; k++)
//                    gta_popP[i-1][k] = gta_popP[j][k];
//                  gta_popC[i-1] = gta_popC[j];
//                  for (k = 0; k < i_D; k++)
//                    gta_popP[j][k] = tempP[k];
//                  gta_popC[j] = tempC;
//                    done = 0;
//                    /* if a swap has been made we are not finished yet */
//              }  /* if */
//          }  /* for */
//        } while (!done);   /* while */
//      } /*while (step > 1) */
//      /* now the best NP are in first NP places in gta_pop, use them */
//      for (i = 0; i < i_NP; i++) {
//        for (j = 0; j < i_D; j++)
//          gta_newP[i][j] = gta_popP[i][j];
//        gta_newC[i] = gta_popC[i];
//      }
    } /*i_bs_flag*/

    /* have selected NP mutants move on to next generation */
    for (i = 0; i < i_NP; i++) {
      for (j = 0; j < i_D; j++)
        ngta_oldP[i+i_NP*j] = ngta_newP[i+i_NP*j];
      ngta_oldC[i] = ngta_newC[i];
    }
    for (j = 0; j < i_D; j++)
	t_bestitP[j] = gt_bestP[j];
      
    if( trace > 0 ) {
      if( (i_iter % trace) == 0 ) {
        Rprintf("Iteration: %d bestvalit: %f bestmemit:", i_iter, t_bestC);
        for (j = 0; j < i_D; j++)
          Rprintf("%12.6f", gt_bestP[j]);
        Rprintf("\n");
      }
    }

    /* check for user interrupt */
    /*if( i_iter % 10000 == 999 ) R_CheckUserInterrupt();*/

    /* check relative tolerance (as in src/main/optim.c) */


    if( gd_bestvalit[i_iter-1] -  t_bestC <
        (d_reltol * (fabs(gd_bestvalit[i_iter-1]) + d_reltol))) {
	i_iter_tol++;
    } else {
      i_iter_tol = 0;
    }

  } /* end iteration loop */

  /* last population */
  k = 0;
  for (i = 0; i < i_NP; i++) {
    for (j = 0; j < i_D; j++) {
      gd_pop[k] = ngta_oldP[i+i_NP*j];
      k++;
    }
  }

  *gi_iter = i_iter;
  *gt_bestC = t_bestC;

  PutRNGstate();
  UNPROTECT(P+1);

}

void permute(int ia_urn2[], int i_urn2_depth, int i_NP, int i_avoid, int ia_urn1[])
/********************************************************************
 ** Function       : void permute(int ia_urn2[], int i_urn2_depth)
 ** Author         : Rainer Storn (w/bug fixes contributed by DEoptim users)
 ** Description    : Generates i_urn2_depth random indices ex [0, i_NP-1]
 **                  which are all distinct. This is done by using a
 **                  permutation algorithm called the "urn algorithm"
 **                  which goes back to C.L.Robinson.
 ** Functions      : -
 ** Globals        : -
 ** Parameters     : ia_urn2       (O)    array containing the random indices
 **                  i_urn2_depth  (I)    number of random indices (avoided index included)
 **                  i_NP          (I)    range of indices is [0, i_NP-1]
 **                  i_avoid       (I)    is the index to avoid and is located in
 **                                       ia_urn2[0].
 ** Preconditions  : # Make sure that ia_urn2[] has a length of i_urn2_depth.
 **                  # i_urn2_depth must be smaller than i_NP.
 ** Postconditions : # the index to be avoided is in ia_urn2[0], so fetch the
 **                   indices from ia_urn2[i], i = 1, 2, 3, ..., i_urn2_depth.
 ** Return Value   : -
 *********************************************************************/
{
  GetRNGstate();

  int k = i_NP;
  int i_urn1 = 0;
  int i_urn2 = 0;
  for (int i = 0; i < i_NP; i++)
    ia_urn1[i] = i; /* initialize urn1 */

  i_urn1 = i_avoid;                    /* get rid of the index to be avoided and place it in position 0. */
  while (k > i_NP - i_urn2_depth)      /* i_urn2_depth is the amount of indices wanted (must be <= NP) */
  {
    ia_urn2[i_urn2] = ia_urn1[i_urn1]; /* move it into urn2 */
    ia_urn1[i_urn1] = ia_urn1[k-1];    /* move highest index to fill gap */
    k = k - 1;                         /* reduce number of accessible indices */
    i_urn2 = i_urn2 + 1;               /* next position in urn2 */
    i_urn1 = (int)(unif_rand() * k);   /* choose a random index */
  }

  PutRNGstate();

}
