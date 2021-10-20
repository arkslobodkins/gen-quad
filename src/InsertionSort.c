/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */

#include "InsertionSort.h"
#include <stdlib.h>

// Sorts current entries of array in ascending order
void InsertionSort(int num_entries, double *array, int *arrayIndex)
{

   double temp = -1.0;
   // Sort s in ascending order
   for(int i = 0; i < num_entries; ++i)
   {
      temp = array[i];
      int j = i;
      while( (j > 0) && (array[j-1] > temp) )
      {
         array[j] = array[j-1];
         arrayIndex[j] = arrayIndex[j-1];
         --j;
      }
      array[j] = temp;
      arrayIndex[j] = i;
   }
}


