//------Prevent multiple includes of de.h----------------------------------
#ifndef _DE_H
#define _DE_H

//------General constants--------------------------------------------------
#define MAXDIM   20        // maximum number of dimensions i.e. parameters. 
#define MAXPOP  2000       // number of random vectors to be stored. Watch  
			               // out! gi_D must be <= 33 because w_index writes 
			               // up to location 3*gi_D-1.    
#define MAXCOST  1    // maximum number of objectives to be minimized
#define MAXCONST 20   // maximum number of constraints

#define TRUE  1
#define FALSE 0

//------Typedefs---------------------------------------------------
typedef struct t_pop 
//*************************************
//** Definition of population member
//*************************************
{
   float fa_vector[MAXDIM];         //parameter vector
   float fa_cost[MAXCOST];          //vector of objectives (costs)
   float fa_constraint[MAXCONST];   //vector of constraints
} t_pop;

#endif // _DE_H
