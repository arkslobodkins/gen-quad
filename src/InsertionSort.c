/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */

#include "InsertionSort.h"
#include <stdlib.h>

// Sorts current indices of an array according to its corresponding norm
void InsertionSort(int num_entries, double *norms, int *arrayIndex)
{
   for(int i = 0; i < num_entries; ++i)
   {
      int j = i;
      double temp = norms[i];
      while( (j > 0) && (norms[j-1] > temp) )
      {
         arrayIndex[j] = arrayIndex[j-1];
         norms[j] = norms[j-1];
         --j;
      }
      arrayIndex[j] = i;
      norms[j] = temp;
   }
}


