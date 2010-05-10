#include "de.h"


/**C*F****************************************************************
**                                                                  
** Function       :void sort (t_pop ta_ary[], int i_len)                                        
**                                                                  
** Author         :Rainer Storn                                     
**                                                                  
** Description    :Shell-sort procedure which sorts array ta_ary[] according
**                 to ta_ary[].fa_cost[0] in ascending order.                 
**                                                                  
** Functions      :-                                                
**                                                                  
** Globals        :-
**                                                                
**                                                                  
** Parameters     :ta_ary[]    (I/O)   population array 
**                 i_len        (I)    length of array to be sorteds   
**                                                                  
** Preconditions  :-                     
**                                                                  
** Postconditions :ta_ary[] will be sorted in ascending order (according to fa_cost[0]) 
**
** Return Value   :-                                            
**                                                                  
***C*F*E*************************************************************/
void sort (t_pop ta_ary[], int i_len)
{
  int   done;
  int   step, bound, i, j;
  t_pop temp;

  step = i_len;  //array length
  while (step > 1) 
  {
    step /= 2;	//halve the step size
    do 
    {
      done   = TRUE;
      bound  = i_len - step;
      for (j = 0; j < bound; j++) 
      {
	    i = j + step + 1;
	    if (ta_ary[j].fa_cost[0] > ta_ary[i-1].fa_cost[0]) 	
	    {
	       temp     = ta_ary[i-1];
	       ta_ary[i-1] = ta_ary[j];
	       ta_ary[j]   = temp;
	       done = FALSE; //if a swap has been made we are not finished yet
	    }  // if
      }  // for
    } while (done == FALSE);   // while
  } //while (step > 1)
} //end of sort()
