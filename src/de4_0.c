/***************************************************************

Implementation of DE based loosely on DE-Engine v4.0, Rainer Storn, 2004   
by Katharine Mullen, 2009

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

//------Include section----------------------------------------------

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rdefines.h>

SEXP getListElement(SEXP list, char *str);
SEXP deoptim(SEXP lower, SEXP upper, SEXP fn, SEXP control, SEXP rho);
void devol(double VTR, double f_weight, double fcross, int i_bs_flag, 
           double *lower, double *upper, SEXP fcall, SEXP rho, int i_trace,
	   int i_specinitialpop, double *initialpopv, int i_storepopfrom, 
	   int i_storepopfreq, int i_check_winner, int i_av_winner);
void permute(int ia_urn2[], int i_urn2_depth, int i_NP, int i_avoid);
double evaluate(long *l_nfeval, double *param, int i_D, SEXP fcall, SEXP env);

//------Global variables--------------------------------------------
// note: these global variables are not necessary - artifact of Storn DE implement
// that should be removed in future 

double VTR;               // value to reach
int   gi_iter;            // generation counter
int   gi_strategy;        // chooses DE-strategy
long  gl_nfeval;          // number of function evaluations      
int   gi_D;               // Dimension of parameter vector
int   gi_NP;              // Number of population members
int   gi_itermax;         // Maximum number of generations

// Data structures for parameter vectors 
double **gta_popP;  
double **gta_oldP;  
double **gta_newP;  
double **gta_swapP; 
double *gt_tmpP;
double *gt_bestP; 

// Data structures for objective function values 
// associated with parameter vectors
double *gta_popC;  
double *gta_oldC;  
double *gta_newC;  
double *gta_swapC; 
double gt_tmpC;
double gt_bestC; 
      
double *gd_pop, *gd_storepop, *gd_bestmemit, *gd_bestvalit;

//------General functions-----------------------------------------

SEXP DEoptimC(SEXP lower, SEXP upper, SEXP fn, SEXP control, SEXP rho)
{
  int i, j;
  int i_bs_flag, i_trace, i_storepopfrom, i_storepopfreq, i_specinitialpop, 
    i_check_winner, i_av_winner;
  double f_cross, f_weight;
  SEXP sexp_bestmem, sexp_bestval, sexp_nfeval, sexp_iter,
    out, out_names, sexp_pop, sexp_storepop, sexp_bestmemit, sexp_bestvalit;

  if (!isFunction(fn))
    error("fn is not a function!");
  if (!isEnvironment(rho))
    error("rho is not an environment!");

  //-----Initialization of annealing parameters-------------------------
  VTR = (double)NUMERIC_VALUE(getListElement(control, "VTR"));
  gi_strategy = (int)NUMERIC_VALUE(getListElement(control, "strategy"));
  gi_itermax = (int)NUMERIC_VALUE(getListElement(control, "itermax"));
  gl_nfeval = (long)NUMERIC_VALUE(getListElement(control, "nfeval"));
  gi_D = (int)NUMERIC_VALUE(getListElement(control, "npar"));
  gi_NP = (int)NUMERIC_VALUE(getListElement(control, "NP"));
  i_storepopfrom = (int)NUMERIC_VALUE(getListElement(control, "storepopfrom"))-1;
  i_storepopfreq = (int)NUMERIC_VALUE(getListElement(control, "storepopfreq"));
  i_specinitialpop = (int)NUMERIC_VALUE(getListElement(control, "specinitialpop"));

  double *initialpopv = (double *)NUMERIC_POINTER(getListElement(control, "initialpop"));

  double f_lower[gi_D];
  double f_upper[gi_D];

  double *lowd = NUMERIC_POINTER(lower);
  double *uppd = NUMERIC_POINTER(upper);

  for (i = 0; i < gi_D; i++) {
    f_lower[i] = (double)(lowd[i]);
    f_upper[i] = (double)(uppd[i]);
  }

  f_weight = NUMERIC_VALUE(getListElement(control, "F"));
  f_cross = NUMERIC_VALUE(getListElement(control, "CR"));
  i_bs_flag = NUMERIC_VALUE(getListElement(control, "bs"));
  i_trace = NUMERIC_VALUE(getListElement(control, "trace"));
  i_check_winner = NUMERIC_VALUE(getListElement(control, "checkWinner"));
  i_av_winner = NUMERIC_VALUE(getListElement(control, "avWinner"));

  //---optimization--------------------------------------
  devol(VTR, f_weight, f_cross, i_bs_flag, f_lower, f_upper, fn, rho, i_trace,
	i_specinitialpop, initialpopv, i_storepopfrom, i_storepopfreq, 
	i_check_winner, i_av_winner);
  //---end optimization----------------------------------

  PROTECT(sexp_bestmem = NEW_NUMERIC(gi_D));
  for (i = 0; i < gi_D; i++) {
    NUMERIC_POINTER(sexp_bestmem)[i] = gt_bestP[i];
  }

  j = gi_NP * gi_D;
  PROTECT(sexp_pop = NEW_NUMERIC(j));
  for (i = 0; i < j; i++)
    NUMERIC_POINTER(sexp_pop)[i] = gd_pop[i];
      
  j =  ceil((gi_itermax - i_storepopfrom) / i_storepopfreq) * gi_NP * gi_D;
  PROTECT(sexp_storepop = NEW_NUMERIC(j));
  for (i = 0; i < j; i++)
    NUMERIC_POINTER(sexp_storepop)[i] = gd_storepop[i];

  j = gi_iter * gi_D;
  PROTECT(sexp_bestmemit = NEW_NUMERIC(j));
  for (i = 0; i < j; i++) 
    NUMERIC_POINTER(sexp_bestmemit)[i] = gd_bestmemit[i];
  j = gi_iter;
  PROTECT(sexp_bestvalit = NEW_NUMERIC(j));
  for (i = 0; i < j; i++)
    NUMERIC_POINTER(sexp_bestvalit)[i] = gd_bestvalit[i];

  PROTECT(sexp_bestval = NEW_NUMERIC(1));
  NUMERIC_POINTER(sexp_bestval)[0] = gt_bestC;

  PROTECT(sexp_nfeval = NEW_INTEGER(1));
  INTEGER_POINTER(sexp_nfeval)[0] = (int)gl_nfeval;

  PROTECT(sexp_iter = NEW_INTEGER(1));
  INTEGER_POINTER(sexp_iter)[0] = gi_iter;

  PROTECT(out = NEW_LIST(8));
  SET_VECTOR_ELT(out, 0, sexp_bestmem);
  SET_VECTOR_ELT(out, 1, sexp_bestval);
  SET_VECTOR_ELT(out, 2, sexp_nfeval);
  SET_VECTOR_ELT(out, 3, sexp_iter);
  SET_VECTOR_ELT(out, 4, sexp_bestmemit);
  SET_VECTOR_ELT(out, 5, sexp_bestvalit);
  SET_VECTOR_ELT(out, 6, sexp_pop);
  SET_VECTOR_ELT(out, 7, sexp_storepop);

  PROTECT(out_names = NEW_STRING(8));
  SET_STRING_ELT(out_names, 0, mkChar("bestmem"));
  SET_STRING_ELT(out_names, 1, mkChar("bestval"));
  SET_STRING_ELT(out_names, 2, mkChar("nfeval"));
  SET_STRING_ELT(out_names, 3, mkChar("iter"));
  SET_STRING_ELT(out_names, 4, mkChar("bestmemit"));
  SET_STRING_ELT(out_names, 5, mkChar("bestvalit"));
  SET_STRING_ELT(out_names, 6, mkChar("pop"));
  SET_STRING_ELT(out_names, 7, mkChar("storepop"));

  SET_NAMES(out, out_names);

  UNPROTECT(10);

  return out;
}

void devol(double VTR, double f_weight, double f_cross, int i_bs_flag, double *lower,
	   double *upper, SEXP fcall, SEXP rho, int trace, int i_specinitialpop,
	   double *initialpopv, int i_storepopfrom, int i_storepopfreq, 
	   int i_check_winner, int i_av_winner)
{
#define URN_DEPTH  5   //4 + one index to avoid

  int i, j, k, x;  // counting variables
  int i_r1, i_r2, i_r3, i_r4;  // placeholders for random indexes

  int i_itermax = gi_itermax;
  int ia_urn2[URN_DEPTH];
  int i_nstorepop, i_xav;
  i_nstorepop = ceil((gi_itermax - i_storepopfrom) / i_storepopfreq);
  
  int popcnt, bestacnt, same; //lazy cnters

  double *fa_minbound = lower;
  double *fa_maxbound = upper;
  double f_jitter, f_dither;
 
  double *t_bestitP;
  double *t_tmpP; 

  double t_bestitC;
  double t_tmpC, tmp_best; 
  
  double p[gi_NP][gi_D];
  double m[gi_D];
  
  // vars for sort 
  int i_len;
  int   done;
  int   step, bound;
  double *tempP;
  double tempC;

  double besta[gi_itermax * gi_D];
  double bestva[gi_itermax];
  double storepop[i_nstorepop * gi_NP * gi_D];
  double pop[gi_NP * gi_D];
  double initialpop[gi_NP][gi_D];

  gta_popP = (double **)R_alloc(gi_NP*2,sizeof(double *));
  for (int i = 0; i < (gi_NP*2); i++) 
    gta_popP[i] = (double *)R_alloc(gi_D,sizeof(double));
  gta_oldP = (double **)R_alloc(gi_NP,sizeof(double *));
  for (int i = 0; i < gi_NP; i++) 
    gta_oldP[i] = (double *)R_alloc(gi_D,sizeof(double));
  gta_newP = (double **)R_alloc(gi_NP,sizeof(double *));
  for (int i = 0; i < gi_NP; i++) 
    gta_newP[i] = (double *)R_alloc(gi_D,sizeof(double));
  
  gta_popC = (double *)R_alloc(gi_NP*2,sizeof(double));
  gta_oldC = (double *)R_alloc(gi_NP,sizeof(double));
  gta_newC = (double *)R_alloc(gi_NP,sizeof(double));
  
  gt_bestP = (double *)R_alloc(1,sizeof(double) * gi_D);
  gt_tmpP = (double *)R_alloc(1,sizeof(double) * gi_D);
  
  t_bestitP = (double *)R_alloc(1,sizeof(double) * gi_D);
  t_tmpP = (double *)R_alloc(1,sizeof(double) * gi_D);
  tempP = (double *)R_alloc(1,sizeof(double) * gi_D);
  
  GetRNGstate();

  gta_popP[0][0] = 0;
 
  // initialize initial popuplation
  for (int i = 0; i < gi_NP; i++) {
    for (int j = 0; j < gi_D; j++) {
      initialpop[i][j] = 0.0;
    }
  }

  // initialize best members
  for (int i = 0; i <= gi_itermax * gi_D; i++)
    besta[i] = 0.0;

  // initialize best values
  for (int i = 0; i <= gi_itermax; i++)
    bestva[i] = 0.0;

  // initialize best population
  for (int i = 0; i < gi_NP * gi_D; i++)
    pop[i] = 0.0;

  // initialize stored populations
  if (i_nstorepop < 0)
    i_nstorepop = 0;

  for (int i = 0; i < (i_nstorepop * gi_NP * gi_D); i++)
    storepop[i] = 0.0;
      
  // if initial population provided, initialize with values
  if (i_specinitialpop > 0) {
    k = 0;
    
    for (j = 0; j < gi_D; j++) {
      for (i = 0; i < gi_NP; i++) {
	initialpop[i][j] = initialpopv[k];
	k += 1;
      }
    }
  }
  gl_nfeval = 0;  // number of function evaluations
  
  //------Initialization-----------------------------
  
  if (i_specinitialpop <= 0) { // random initial member
    for (j = 0; j < gi_D; j++) 
      gta_popP[0][j] = fa_minbound[j] +
	unif_rand() * (fa_maxbound[j] - fa_minbound[j]);
  }
  else // or user-specified initial member
    for (j = 0; j < gi_D; j++)
      gta_popP[0][j] = initialpop[0][j];
 
  gta_popC[0] = evaluate(&gl_nfeval, gta_popP[0], gi_D, fcall, rho); 
  
  gt_bestC=gta_popC[0];
  for (j = 0; j < gi_D; j++)  
    gt_bestP[j]=gta_popP[0][j];

  for (i = 1; i < gi_NP; i++) {
    for (j = 0; j < gi_D; j++) {
      if (i_specinitialpop <= 0) {
	gta_popP[i][j] = fa_minbound[j] +
	  unif_rand() * (fa_maxbound[j] - fa_minbound[j]);

      }
      else
	gta_popP[i][j] = initialpop[i][j];
    } 
    gta_popC[i] = evaluate(&gl_nfeval, gta_popP[i], gi_D, fcall, rho);
   
    if (gta_popC[i] <= gt_bestC) {
      gt_bestC = gta_popC[i];
      for (j = 0; j < gi_D; j++)  
	gt_bestP[j]=gta_popP[i][j];
    }
  }
 
  
  //---assign pointers to current ("old") population---
  gta_oldP = gta_popP;
  gta_oldC = gta_popC;
  
  //------Iteration loop--------------------------------------------
  gi_iter = 0;
  popcnt = 0;
  bestacnt = 0;
  i_xav = 1;
  
  // loop
  while ((gi_iter < i_itermax) && (gt_bestC > VTR))
    {
      // store intermediate populations
      if (gi_iter % i_storepopfreq == 0 && gi_iter >= i_storepopfrom) {
	for (i = 0; i < gi_NP; i++) {
	  for (j = 0; j < gi_D; j++) {
	    storepop[popcnt] = gta_oldP[i][j];
	    popcnt++;
	  }
	}
       } // end store pop
      
      //store the best member
      for(j = 0; j < gi_D; j++) {
	besta[bestacnt] = gt_bestP[j];
	bestacnt++;
      }
      // store the best value
      bestva[gi_iter] = gt_bestC;
      
      for (j = 0; j < gi_D; j++) 
      	t_bestitP[j] = gt_bestP[j];
      //t_bestitP = gt_bestP;
      t_bestitC = gt_bestC;
      
      gi_iter++;
     
      //----computer dithering factor -----------------
      f_dither = f_weight + unif_rand() * (1.0 - f_weight);

      //----start of loop through ensemble------------------------
      for (i = 0; i < gi_NP; i++) {

	//t_tmpP is the vector to mutate and eventually select
	for (j = 0; j < gi_D; j++) 
	  t_tmpP[j] = gta_oldP[i][j];
	t_tmpC = gta_oldC[i];

	permute(ia_urn2, URN_DEPTH, gi_NP, i); //Pick 4 random and distinct

	i_r1 = ia_urn2[1];  //population members
	i_r2 = ia_urn2[2];
	i_r3 = ia_urn2[3];
	i_r4 = ia_urn2[4];
	
	//===Choice of strategy=======================================================
	//---classical strategy DE/rand/1/bin-----------------------------------------
	if (gi_strategy == 1) {
	  
	  j = (int)(unif_rand() * gi_D); // random parameter
	  k = 0;
	  do {
	    // add fluctuation to random target
	    t_tmpP[j] = gta_oldP[i_r1][j] +
	      f_weight * (gta_oldP[i_r2][j] - gta_oldP[i_r3][j]);

	    j = (j + 1) % gi_D;
	    k++;
	  }while((unif_rand() < f_cross) && (k < gi_D));

	}
	//---DE/local-to-best/1/bin---------------------------------------------------
	else if (gi_strategy == 2) {
	 
	  j = (int)(unif_rand() * gi_D); // random parameter
	  k = 0;
	  do {
	    // add fluctuation to random target
	   
	    t_tmpP[j] = t_tmpP[j] + 
	      f_weight * (t_bestitP[j] - t_tmpP[j]) +
	      f_weight * (gta_oldP[i_r2][j] - gta_oldP[i_r3][j]);
	    j = (j + 1) % gi_D;
	    k++;
	  }while((unif_rand() < f_cross) && (k < gi_D));

	}
	//---DE/best/1/bin with jitter------------------------------------------------
	else if (gi_strategy == 3) {
	 	  
	  j = (int)(unif_rand() * gi_D); // random parameter
	  k = 0;
	  do {
	    // add fluctuation to random target
	    f_jitter = 0.0001 * unif_rand() + f_weight;
	    t_tmpP[j] = t_bestitP[j] +
	      f_jitter * (gta_oldP[i_r1][j] - gta_oldP[i_r2][j]);
	    
	    j = (j + 1) % gi_D;
	    k++;
	  }while((unif_rand() < f_cross) && (k < gi_D));

	}
	//---DE/rand/1/bin with per-vector-dither-------------------------------------
	else if (gi_strategy == 4) {
		  
	  j = (int)(unif_rand() * gi_D); // random parameter
	  k = 0;
	  do {
	    // add fluctuation to random target
	    t_tmpP[j] = gta_oldP[i_r1][j] +
	      (f_weight + unif_rand()*(1.0 - f_weight))*
	      (gta_oldP[i_r2][j]-gta_oldP[i_r3][j]);

	    j = (j + 1) % gi_D;
	    k++;
	  }while((unif_rand() < f_cross) && (k < gi_D));

	}
	//---DE/rand/1/bin with per-generation-dither---------------------------------
	else if (gi_strategy == 5) {
	  
	  j = (int)(unif_rand() * gi_D); // random parameter
	  k = 0;
	  do {
	    // add fluctuation to random target
	    t_tmpP[j] = gta_oldP[i_r1][j] +
	      f_dither * (gta_oldP[i_r2][j] - gta_oldP[i_r3][j]);
	    
	    j = (j + 1) % gi_D;
	    k++;
	  }while((unif_rand() < f_cross) && (k < gi_D));
       
	}
	//---variation to DE/rand/1/bin: either-or-algorithm--------------------------
	else {
	  
	  j = (int)(unif_rand() * gi_D); // random parameter
	  k = 0;
	  if (unif_rand() < 0.5) { //differential mutation, Pmu = 0.5
	    do {
	      // add fluctuation to random target
	      t_tmpP[j] = gta_oldP[i_r1][j] +
		f_weight * (gta_oldP[i_r2][j] - gta_oldP[i_r3][j]);
	      
	      j = (j + 1) % gi_D;
	      k++;
	    }while((unif_rand() < f_cross) && (k < gi_D));
	  }
	  else {
	    //recombination with K = 0.5*(F+1) -. F-K-Rule
	    do {
	      // add fluctuation to random target
	      t_tmpP[j] = gta_oldP[i_r1][j] +
		0.5 * (f_weight + 1.0) * (gta_oldP[i_r2][j]
					  + gta_oldP[i_r3][j] - 2 * gta_oldP[i_r1][j]);
	      
	      j = (j + 1) % gi_D;
	      k++;
	    }while((unif_rand() < f_cross) && (k < gi_D));

	  }
	}//end if (gi_strategy ...
	
	//----boundary constraints, bounce-back method was not enforcing bounds correctly
	for (j = 0; j < gi_D; j++) {
	  if (t_tmpP[j] < fa_minbound[j]) {
	    t_tmpP[j] = fa_minbound[j] +
	      unif_rand() * (fa_maxbound[j] - fa_minbound[j]);
	  }
	  if (t_tmpP[j] > fa_maxbound[j]) {
	    t_tmpP[j] =  fa_maxbound[j] -
	      unif_rand() * (fa_maxbound[j] - fa_minbound[j]);
	  }
	}

	//------Trial mutation now in t_tmpP-----------------
	// Evaluate mutant in t_tmpP[]

	t_tmpC = evaluate(&gl_nfeval, t_tmpP, gi_D, fcall, rho); 
	
	// note that i_bs_flag means that we will choose the 
	//best NP vectors from the old and new population later 
	if (t_tmpC <= gta_oldC[i] || i_bs_flag) {
	  // replace target with mutant
	  for (j = 0; j < gi_D; j++) 
	    gta_newP[i][j]=t_tmpP[j];
	  gta_newC[i]=t_tmpC;
	  if (t_tmpC <= gt_bestC) {
	    for (j = 0; j < gi_D; j++) 
	      gt_bestP[j]=t_tmpP[j];
	    gt_bestC=t_tmpC;
	  }
	} 
	else {
	  for (j = 0; j < gi_D; j++) 
	    gta_newP[i][j]=gta_oldP[i][j];
	  gta_newC[i]=gta_oldC[i];
	  
	}
      } // End mutation loop through pop.
     
     
      if(i_bs_flag) { 
	// examine old and new pop. and take the best NP members
	// into next generation 
	for (i = 0; i < gi_NP; i++) {
	  for (j = 0; j < gi_D; j++) 
	    gta_popP[i][j] = gta_oldP[i][j];
	  gta_popC[i] = gta_oldC[i];
	}
	for (i = 0; i < gi_NP; i++) {
	  for (j = 0; j < gi_D; j++) 
	    gta_popP[gi_NP+i][j] = gta_newP[i][j];
	  gta_popC[gi_NP+i] = gta_newC[i];
	}
	i_len = 2 * gi_NP;
	step = i_len;  //array length
	while (step > 1) {
	  step /= 2;	//halve the step size
	  do {
	    done = TRUE;
	    bound  = i_len - step;
	    for (j = 0; j < bound; j++) 
	      {
		i = j + step + 1;
		if (gta_popC[j] > gta_popC[i-1]) 	
		  {
		    for (k = 0; k < gi_D; k++) 
		      tempP[k] = gta_popP[i-1][k];
		    tempC = gta_popC[i-1];
		    for (k = 0; k < gi_D; k++) 
		      gta_popP[i-1][k] = gta_popP[j][k];
		    gta_popC[i-1] = gta_popC[j];
		    for (k = 0; k < gi_D; k++) 
		      gta_popP[j][k] = tempP[k];
		    gta_popC[j] = tempC;
		      done = FALSE; 
		      //if a swap has been made we are not finished yet
		  }  // if
	      }  // for
	  } while (!done);   // while
	} //while (step > 1)
	// now the best NP are in first NP places in gta_pop, use them 
	for (i = 0; i < gi_NP; i++) {
	  for (j = 0; j < gi_D; j++) 
	    gta_newP[i][j] = gta_popP[i][j];
	  gta_newC[i] = gta_popC[i];
	  
	}
	// make sure best is best of old and new gen.
	if (gta_newC[0] <= gt_bestC) {
	  for (j = 0; j < gi_D; j++) 
	    gt_bestP[j]=gta_newP[i][j];
	  gt_bestC=gta_newC[0];
	}
      }
      // have selected NP mutants
      // move on to next generation
      for (i = 0; i < gi_NP; i++) {
	for (j = 0; j < gi_D; j++) 
	  gta_oldP[i][j] = gta_newP[i][j];
	gta_oldC[i] = gta_newC[i];
      }
      //check if the best stayed the same, if necessary
      if(i_check_winner)  {
	same = 1;
	for (j = 0; j < gi_D; j++)
	  if(t_bestitP[j] != gt_bestP[j]) {
	    same = 0;
	  }
	if(same && gi_iter > 1)  {
	  i_xav++;
	  //if re-evaluation of winner 
	  tmp_best = evaluate(&gl_nfeval, gt_bestP, gi_D, fcall, rho);
	 
	  //possibly letting the winner be the average of all past
	  //generations
	  if(i_av_winner)
	    gt_bestC = ((1/(double)i_xav) * gt_bestC) 
	      + ((1/(double)i_xav) * tmp_best) + (bestva[gi_iter-1] * ((double)(i_xav - 2))/(double)i_xav);
	  else
	    gt_bestC = tmp_best;
	
	}
	else {
	  i_xav = 1;
	}
	
      }
      for (j = 0; j < gi_D; j++) 
	t_bestitP[j] = gt_bestP[j];
      t_bestitC = gt_bestC;
      
      if(trace > 0) {
	Rprintf("Iteration: %d bestvalit: %f bestmemit:", gi_iter, gt_bestC);
	for (j = 0; j < gi_D; j++)
	  Rprintf("%12.6f", gt_bestP[j]);
	Rprintf("\n");
      }
    } // end loop through generations
  
  // last population
  k = 0;
  for (i = 0; i < gi_NP; i++) {
    for (j = 0; j < gi_D; j++) {
      pop[k] = gta_oldP[i][j];      
      k++;
    }
  }

  // global assignment
  gd_bestmemit = besta;
  gd_bestvalit = bestva;
  gd_pop = pop;
  gd_storepop = storepop;
  
  PutRNGstate();
    
}

void permute(int ia_urn2[], int i_urn2_depth, int i_NP, int i_avoid)
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
  int  i, k, i_urn1, i_urn2;
  int  *ia_urn1;      //urn holding all indices
  ia_urn1 = (int *)R_alloc(gi_NP,sizeof(int));
  for(i = 0; i < gi_NP; i++)
    ia_urn1[i] = 0;
  
  GetRNGstate();
	  
  k = i_NP;
  i_urn1 = 0;
  i_urn2 = 0;
  for (i = 0; i < i_NP; i++)
    ia_urn1[i] = i; //initialize urn1

  i_urn1 = i_avoid;                  //get rid of the index to be avoided and place it in position 0.
  while (k > i_NP - i_urn2_depth)     //i_urn2_depth is the amount of indices wanted (must be <= NP)
    {
      ia_urn2[i_urn2] = ia_urn1[i_urn1];      //move it into urn2
      ia_urn1[i_urn1] = ia_urn1[k-1]; //move highest index to fill gap
      k = k - 1;                        //reduce number of accessible indices
      i_urn2 = i_urn2 + 1;            //next position in urn2
      i_urn1 = (int)(unif_rand() * k);    //choose a random index
    }

  PutRNGstate();
	
}
