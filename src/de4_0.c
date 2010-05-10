/***************************************************************

Modification of DE-Engine v4.0, Rainer Storn, 2004   
by Katharine Mullen, 2009

Allows compilation in an R package

Storn's MS Visual C++ v5.0 was downloaded from 
http://www.icsi.berkeley.edu/~storn/DeWin.zip
and translated to vanilla C, modified to input and output from R,
store the population at each generation, etc.
***************************************************************/

//------Include section----------------------------------------------

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "de.h"
#include <R.h>
#include <Rdefines.h>

SEXP getListElement(SEXP list, char *str);
SEXP deoptim(SEXP lower, SEXP upper, SEXP fn, SEXP control, SEXP rho);
void devol(float VTR, float f_weight, float fcross, int i_bs_flag, 
           float *lower, float *upper, SEXP fcall, SEXP rho, int i_trace,
	   int i_specinitialpop, double *initialpopv, int i_storepopfrom, 
	   int i_storepopfreq);
void sort (t_pop ary[], int len);
void assigna2b(int D, float a[], float b[]);
void permute(int ia_urn2[], int i_urn2_depth, int i_NP, int i_avoid);
t_pop evaluate(t_pop t_tmp, long *l_nfeval, t_pop *tpa_array, int i_NP,
               SEXP fcall, SEXP env);
int left_vector_wins(t_pop t_trial, t_pop t_target);

//------Global variables--------------------------------------------

float VTR;                 // value to reach
int   gi_iter;             // generation counter
int   gi_strategy;        // chooses DE-strategy
long  gl_nfeval;          // number of function evaluations      
int   gi_D;               // Dimension of parameter vector
int   gi_NP;              // Number of population members
int   gi_itermax;         // Maximum number of generations
t_pop gta_pop[2*MAXPOP];  // the two populations are put into one array side by side.      
t_pop gt_best;            // current best population member
t_pop *gpta_old, *gpta_new, *gpta_swap;
double *gd_pop, *gd_storepop, *gd_bestmemit, *gd_bestvalit;
int gi_nstorepop;

//------General functions-----------------------------------------

SEXP DEoptimC(SEXP lower, SEXP upper, SEXP fn, SEXP control, SEXP rho)
{
   int i, j;
   int i_bs_flag, i_trace, i_storepopfrom, i_storepopfreq, i_specinitialpop;
   float f_cross, f_weight;
   SEXP sexp_bestmem, sexp_bestval, sexp_nfeval, sexp_iter,
        out, out_names, sexp_pop, sexp_storepop, sexp_bestmemit, sexp_bestvalit;

   if (!isFunction(fn))
      error("fn is not a function!");
   if (!isEnvironment(rho))
      error("rho is not an environment!");

   //-----Initialization of annealing parameters-------------------------
   VTR = (float)NUMERIC_VALUE(getListElement(control, "VTR"));
   gi_strategy = (int)NUMERIC_VALUE(getListElement(control, "strategy"));
   gi_itermax = (int)NUMERIC_VALUE(getListElement(control, "itermax"));
   gl_nfeval = (long)NUMERIC_VALUE(getListElement(control, "nfeval"));
   gi_D = (int)NUMERIC_VALUE(getListElement(control, "npar"));
   gi_NP = (int)NUMERIC_VALUE(getListElement(control, "NP"));
   i_storepopfrom = (int)NUMERIC_VALUE(getListElement(control, "storepopfrom")) - 1;
   i_storepopfreq = (int)NUMERIC_VALUE(getListElement(control, "storepopfreq"));
   i_specinitialpop = (int)NUMERIC_VALUE(getListElement(control, "specinitialpop"));
   
   double *initialpopv = (double *)NUMERIC_POINTER(getListElement(control, "initialpop"));

   float f_lower[gi_D];
   float f_upper[gi_D];

   double *lowd = NUMERIC_POINTER(lower);
   double *uppd = NUMERIC_POINTER(upper);

   for (i = 0; i < gi_D; i++) {
      f_lower[i] = (float)(lowd[i]);
      f_upper[i] = (float)(uppd[i]);
   }

   if (gi_D > MAXDIM)
      error("too many parameters!");
   if (gi_NP > MAXPOP)
      error("too many points!");

   f_weight = NUMERIC_VALUE(getListElement(control, "F"));
   f_cross = NUMERIC_VALUE(getListElement(control, "CR"));
   i_bs_flag = NUMERIC_VALUE(getListElement(control, "i_bs_flag"));
   i_trace = NUMERIC_VALUE(getListElement(control, "trace"));

   //---optimization--------------------------------------
   devol(VTR, f_weight, f_cross, i_bs_flag, f_lower, f_upper, fn, rho, i_trace,
	      i_specinitialpop, initialpopv, i_storepopfrom, i_storepopfreq);
   //---end optimization----------------------------------

   PROTECT(sexp_bestmem = NEW_NUMERIC(gi_D));
   for (i = 0; i < gi_D; i++) {
      NUMERIC_POINTER(sexp_bestmem)[i] = gt_best.fa_vector[i];
   }

   j = gi_NP * gi_D;
   PROTECT(sexp_pop = NEW_NUMERIC(j));
   for (i = 0; i < j; i++)
      NUMERIC_POINTER(sexp_pop)[i] = gd_pop[i];
      
   j = gi_nstorepop * gi_NP * gi_D;
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
   NUMERIC_POINTER(sexp_bestval)[0] = gt_best.fa_cost[0];

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

void  assigna2b(int i_D, float fa_a[], float fa_b[])
{
   int j;
   for (j = 0; j < i_D; j++)
      fa_b[j] = fa_a[j];
}


void devol(float VTR, float f_weight, float f_cross, int i_bs_flag, float *lower,
	        float *upper, SEXP fcall, SEXP rho, int trace, int i_specinitialpop,
	        double *initialpopv, int i_storepopfrom, int i_storepopfreq)
{
#define URN_DEPTH  5   //4 + one index to avoid
   //-----Variable declarations--------------------------------------- --
   int i, j, k;  // counting variables
   int i_r1, i_r2, i_r3, i_r4;  // placeholders for random indexes

   int i_itermax = gi_itermax;
   int ia_urn2[URN_DEPTH];

   int popcnt, bestacnt; //lazy cnters

   float *fa_minbound = lower;
   float *fa_maxbound = upper;
   t_pop t_tmp, t_bestit;

   t_pop t_origin;

   float f_jitter, f_dither;

   GetRNGstate();

   // initialize initial popuplation
   double initialpop[gi_NP][gi_D];
   for (int i = 0; i < gi_NP; i++) {
		 for (int j = 0; j < gi_D; j++) {
		 	  initialpop[i][j] = 0.0;
	     }
   }

   // initialize best members
   double besta[gi_itermax * gi_D];
   for (int i = 0; i <= gi_itermax * gi_D; i++)
      besta[i] = 0.0;

   // initialize best values
   double bestva[gi_itermax];
   for (int i = 0; i <= gi_itermax; i++)
      bestva[i] = 0.0;

   // initialize best population
   double pop[gi_NP * gi_D];
   for (int i = 0; i < gi_NP * gi_D; i++)
      pop[i] = 0.0;

   // initialize stored populations
   gi_nstorepop = ceil((gi_itermax - i_storepopfrom) / i_storepopfreq);

	if (gi_nstorepop < 0)
	   gi_nstorepop = 0;

   double storepop[gi_nstorepop * gi_NP * gi_D];
   for (int i = 0; i < gi_nstorepop * gi_NP * gi_D; i++)
      storepop[i] = 0.0;
   
   
   // if initial population provided, initialize with values
   if (i_specinitialpop > 0) {
      k = 0;
      for (i = 0; i < gi_NP; i++) {
         for (j = 0; j < gi_D; j++) {
	         initialpop[i][j] = initialpopv[k];
	         k += 1;
         }
      }
   }
   gl_nfeval = 0;  // reset number of function evaluations

   //------Initialization-----------------------------
   if (i_specinitialpop <= 0) {
     for (j = 0; j < gi_D; j++)
         gta_pop[0].fa_vector[j] = fa_minbound[j] +
			   unif_rand() * (fa_maxbound[j] - fa_minbound[j]);
   }
   else
      for (j = 0; j < gi_D; j++)
         gta_pop[0].fa_vector[j] = initialpop[0][j];

   gta_pop[0] = evaluate(gta_pop[0], &gl_nfeval, &gta_pop[0], gi_NP, fcall, rho);

   gt_best = gta_pop[0];

   for (i = 1; i < gi_NP; i++) {
      for (j = 0; j < gi_D; j++) {
         if (i_specinitialpop <= 0)
            gta_pop[i].fa_vector[j] = fa_minbound[j] +
				   unif_rand() * (fa_maxbound[j] - fa_minbound[j]);
         else
	         gta_pop[i].fa_vector[j] = initialpop[i][j];
      }
      gta_pop[i] = evaluate(gta_pop[i], &gl_nfeval, &gta_pop[0], gi_NP, fcall, rho);
      if (left_vector_wins(gta_pop[i], gt_best) == TRUE)
         gt_best = gta_pop[i];
   }

   t_bestit = gt_best;
    
   //---assign pointers to current ("old") and new population---
   gpta_old = &gta_pop[0];
   gpta_new = &gta_pop[gi_NP];

   //------Iteration loop--------------------------------------------
   gi_iter = 0;
   popcnt = 0;
   bestacnt = 0;
   // loop
   while ((gi_iter < i_itermax) && (gt_best.fa_cost[0] > VTR))
   {
      // store intermediate populations
      if (gi_iter % i_storepopfreq == 0 && gi_iter >= i_storepopfrom) {
         for (i = 0; i < gi_NP; i++) {
	         for (j = 0; j < gi_D; j++) {
	            storepop[popcnt] = (double)gpta_old[i].fa_vector[j];
	            popcnt++;
	         }
         }
      } // end store pop

	   //store the best member
      for(j = 0; j < gi_D; j++) {
		   besta[bestacnt] = gt_best.fa_vector[j];
		   bestacnt++;
      }
      // store the best value
      bestva[gi_iter] = gt_best.fa_cost[0];

      gi_iter++;

      //----computer dithering factor (if needed)-----------------
      f_dither = f_weight + unif_rand() * (1.0 - f_weight);

      //----start of loop through ensemble------------------------
      for (i = 0; i < gi_NP; i++) {

			permute(ia_urn2, URN_DEPTH, gi_NP, i); //Pick 4 random and distinct

			i_r1 = ia_urn2[1];  //population members
	      i_r2 = ia_urn2[2];
	      i_r3 = ia_urn2[3];
	      i_r4 = ia_urn2[4];

         //===Choice of strategy=======================================================
	      //---classical strategy DE/rand/1/bin-----------------------------------------
			if (gi_strategy == 1) {

				assigna2b(gi_D, gpta_old[i].fa_vector, t_tmp.fa_vector);

				j = (int)(unif_rand() * gi_D); // random parameter
	         k = 0;
            do {
 			      // add fluctuation to random target
		         t_tmp.fa_vector[j] = gpta_old[i_r1].fa_vector[j] +
				      f_weight * (gpta_old[i_r2].fa_vector[j] - gpta_old[i_r3].fa_vector[j]);

               j = (j + 1) % gi_D;
		         k++;
		      }while((unif_rand() < f_cross) && (k < gi_D));

            assigna2b(gi_D, gpta_old[i_r1].fa_vector, t_origin.fa_vector);
         }
	      //---DE/local-to-best/1/bin---------------------------------------------------
	      else if (gi_strategy == 2) {

            assigna2b(gi_D, gpta_old[i].fa_vector, t_tmp.fa_vector);

				j = (int)(unif_rand() * gi_D); // random parameter
         	k = 0;
				do {
					// add fluctuation to random target
		    		t_tmp.fa_vector[j] = t_tmp.fa_vector[j] +
				      f_weight * (t_bestit.fa_vector[j] - t_tmp.fa_vector[j]) +
					      f_weight * (gpta_old[i_r2].fa_vector[j] - gpta_old[i_r3].fa_vector[j]);

					j = (j + 1) % gi_D;
		         k++;
				}while((unif_rand() < f_cross) && (k < gi_D));

				assigna2b(gi_D, t_tmp.fa_vector, t_origin.fa_vector);
         }
	      //---DE/best/1/bin with jitter------------------------------------------------
         else if (gi_strategy == 3) {

		      assigna2b(gi_D, gpta_old[i].fa_vector, t_tmp.fa_vector);

			   j = (int)(unif_rand() * gi_D); // random parameter
		      k = 0;
		      do {
			      // add fluctuation to random target
               f_jitter = 0.0001 * unif_rand() + f_weight;
               t_tmp.fa_vector[j] = t_bestit.fa_vector[j] +
			         f_jitter * (gpta_old[i_r1].fa_vector[j] - gpta_old[i_r2].fa_vector[j]);

               j = (j + 1) % gi_D;
               k++;
	         }while((unif_rand() < f_cross) && (k < gi_D));

		      assigna2b(gi_D, t_tmp.fa_vector, t_origin.fa_vector);

         }
		   //---DE/rand/1/bin with per-vector-dither-------------------------------------
         else if (gi_strategy == 4) {

            assigna2b(gi_D, gpta_old[i].fa_vector, t_tmp.fa_vector);

			   j = (int)(unif_rand() * gi_D); // random parameter
            k = 0;
            do {
			      // add fluctuation to random target
		         t_tmp.fa_vector[j] = gpta_old[i_r1].fa_vector[j] +
		            (f_weight + unif_rand()*(1.0 - f_weight))*
		               (gpta_old[i_r2].fa_vector[j]-gpta_old[i_r3].fa_vector[j]);

               j = (j + 1) % gi_D;
		         k++;
		      }while((unif_rand() < f_cross) && (k < gi_D));

		      assigna2b(gi_D, t_tmp.fa_vector, t_origin.fa_vector);
         }
		   //---DE/rand/1/bin with per-generation-dither---------------------------------
	      else if (gi_strategy == 5) {

            assigna2b(gi_D, gpta_old[i].fa_vector, t_tmp.fa_vector);

			   j = (int)(unif_rand() * gi_D); // random parameter
		      k = 0;
            do {
				   // add fluctuation to random target
		         t_tmp.fa_vector[j] = gpta_old[i_r1].fa_vector[j] +
				      f_dither * (gpta_old[i_r2].fa_vector[j] - gpta_old[i_r3].fa_vector[j]);

		         j = (j + 1) % gi_D;
		         k++;
		      }while((unif_rand() < f_cross) && (k < gi_D));

	         assigna2b(gi_D,t_tmp.fa_vector,t_origin.fa_vector);
		   }
		   //---variation to DE/rand/1/bin: either-or-algorithm--------------------------
	      else {
		      assigna2b(gi_D, gpta_old[i].fa_vector, t_tmp.fa_vector);
		      j = (int)(unif_rand() * gi_D); // random parameter
		      k = 0;
		      if (unif_rand() < 0.5) { //differential mutation, Pmu = 0.5
	            do {
				      // add fluctuation to random target
			 	   	t_tmp.fa_vector[j] = gpta_old[i_r1].fa_vector[j] +
				 	      f_weight * (gpta_old[i_r2].fa_vector[j] - gpta_old[i_r3].fa_vector[j]);

			         j = (j + 1) % gi_D;
			         k++;
		         }while((unif_rand() < f_cross) && (k < gi_D));
		      }
		   else {
		      //recombination with K = 0.5*(F+1) --> F-K-Rule
		      do {
			      // add fluctuation to random target
    			   t_tmp.fa_vector[j] = gpta_old[i_r1].fa_vector[j] +
			         0.5 * (f_weight + 1.0) * (gpta_old[i_r2].fa_vector[j]
						   + gpta_old[i_r3].fa_vector[j] - 2 * gpta_old[i_r1].fa_vector[j]);

			      j = (j + 1) % gi_D;
			      k++;
		      }while((unif_rand() < f_cross) && (k < gi_D));
		   }

		      assigna2b(gi_D, gpta_old[i_r1].fa_vector, t_origin.fa_vector);

      }//end if (gi_strategy ...
		
		//----boundary constraints: : removed the 'bounce-back approach' in Storn code
                // to use gpta_old[i].fa_vector[j] instead of the t_origin.fa_vector
			for (j = 0; j < gi_D; j++) {
			  if (t_tmp.fa_vector[j] < fa_minbound[j]) {
			    t_tmp.fa_vector[j] = fa_minbound[j] +
			      unif_rand() * (gpta_old[i].fa_vector[j] - fa_minbound[j]);
              }
			  if (t_tmp.fa_vector[j] > fa_maxbound[j]) {
			    t_tmp.fa_vector[j] = gpta_old[i].fa_vector[j] +
			      unif_rand() * (fa_maxbound[j] - fa_maxbound[j]);
			  }
			}



		//------Trial mutation now in t_tmp-----------------
		t_tmp = evaluate(t_tmp, &gl_nfeval, &gpta_old[0], gi_NP, fcall, rho);  // Evaluate mutant in t_tmp[]

		if (i_bs_flag == TRUE) {
         gpta_new[i] = t_tmp; //save new vector, selection will come later
  		}
		else {
         if (left_vector_wins(t_tmp, gpta_old[i]) == TRUE) {
			   gpta_new[i] = t_tmp; // replace target with mutant

			   if (left_vector_wins(t_tmp, gt_best) == TRUE) {
   			   gt_best = t_tmp;
			   }
		   }
	      else {
			   gpta_new[i] = gpta_old[i];
		   }
      }
	 } // End mutation loop through pop.

    if (i_bs_flag == TRUE) {
	   sort (gpta_old, 2 * gi_NP); //sort array of parents + children
	   gt_best = gpta_old[0];
	 }
    else {
      gpta_swap = gpta_old;
	   gpta_old  = gpta_new;
	   gpta_new  = gpta_swap;
    }

    t_bestit = gt_best;
   if(trace > 0) {
     Rprintf("Iteration: %d bestvalit: %f bestmemit:", gi_iter, gt_best.fa_cost[0]);
      for (i = 0; i < gi_D; i++)
	Rprintf("%12.6f", gt_best.fa_vector[i]);
      Rprintf("\n");
   }
   }

   // last population
   k = 0;
   for (i = 0; i < gi_NP; i++) {
      for (j = 0; j < gi_D; j++) {
         pop[k] = (double)gpta_new[i].fa_vector[j];
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
** Author         : Rainer Storn
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
	int  ia_urn1[MAXPOP] = {0};      //urn holding all indices
	k = i_NP;
	i_urn1 = 0;
	i_urn2 = 0;
	for (i = 0; i < i_NP; i++)
	  ia_urn1[i] = i; //initialize urn1

	i_urn1 = i_avoid;                  //get rid of the index to be avoided and place it in position 0.
	while (k >= i_NP - i_urn2_depth)     //i_urn2_depth is the amount of indices wanted (must be <= NP)
	  {
	   ia_urn2[i_urn2] = ia_urn1[i_urn1];      //move it into urn2
	   ia_urn1[i_urn1] = ia_urn1[k-1]; //move highest index to fill gap
	   k = k - 1;                        //reduce number of accessible indices
	   i_urn2 = i_urn2 + 1;            //next position in urn2
	   i_urn1 = (int)(unif_rand() * k);    //choose a random index
	}
}
