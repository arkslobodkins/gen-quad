/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */

#include "BasisIndices.h"

#include <assert.h>
#include <stdlib.h>

// Recursively computes an array of power multi-indices for
// generating polynomial basis of dimension dim and degree deg.
void BasisIndices(int deg, int dim, INT_8 *f)
{
   assert(dim >= 0 && deg >= 0);

   if(dim == 1) {
      for(int i = 0; i <= deg; ++i) f[i] = i;
      return;
   }

   int counter;
   // compute basis indices using nested recursion if dimension >= 2
   for(int j = 2; counter = 0, j <= dim; ++j)
   {
      for(int i = 0; i <= deg; ++i)
      {
         int size = BasisSize(deg-i, j-1);
         int dimxsize = dim*size;
         INT_8* recursiveF = (INT_8 *)malloc(size*(j-1) * sizeof(INT_8));
         BasisIndices(deg-i, j-1, recursiveF);
         if(j == dim)
         {
            for(int k = 0; k < size; ++k)
            {
               int kxdim = k*dim;
               int kxdim_minus1 = k*(dim-1);
               for(int d = 0; d < dim-1; ++d)
                  f[counter+kxdim+d] = recursiveF[kxdim_minus1+d];
               f[counter+kxdim+j-1] = i;
            }
            counter += dimxsize;
         }
         free(recursiveF);
      }
   }

}


// Computes the minimum number of basis functions to generate
// polynomial basis of dimension dim and degree deg.
int BasisSize(int deg, int dim)
{
   assert(dim >= 0 && deg >= 0);

   unsigned int first = deg + 1;
   unsigned int last = deg + dim;
   unsigned int product = 1;

   unsigned int i = first, counter = 1;
   for(; i <= last; ++i, ++counter)
      product = product * i/counter;

   return product;
}
