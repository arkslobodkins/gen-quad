/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */

#include "BasisIntegrals.h"
#include "BasisIndices.h"
#include "Quadrature.h"

#include <stdlib.h>
#include <string.h>

// Routines that compute analytical integrals of all orthogonal
// basis functions for respective domains and parameters, which
// are all zero due to orthogonality, except for the basis
// function that integrates scalar polynomial "1".

void BasisIntegralsCube(int *dims, int deg, double *integrals)
{
   int dim = dims[0];
   int m = BasisSize(deg, dim);
   memset(integrals, 0, SIZE_DOUBLE(m));
   integrals[0] = 1.0;
}

void BasisIntegralsSimplex(int *dims, int deg, double *integrals)
{
   int dim = dims[0];
   int m = BasisSize(deg, dim);

   memset(integrals, 0, SIZE_DOUBLE(m));
   integrals[0] = 1;
   for(int i = 1; i <= dim; ++i)
      integrals[0] /= i;
}

void BasisIntegralsCubeSimplex(int *dims, int deg, double *integrals)
{
   int dim = dims[0]+dims[1];
   int dim2 = dims[1];
   int m = BasisSize(deg, dim);

   memset(integrals, 0, SIZE_DOUBLE(m));
   integrals[0] = 1;
   for(int i = 1; i <= dim2; ++i)
      integrals[0] /= i;
}

void BasisIntegralsSimplexSimplex(int *dims, int deg, double *integrals)
{
   int dim = dims[0]+dims[1];
   int dim1 = dims[0];
   int dim2 = dims[1];
   int m = BasisSize(deg, dim);

   memset(integrals, 0, SIZE_DOUBLE(m));
   integrals[0] = 1;
   for(int i = 1; i <= dim1; ++i)
      integrals[0] /= i;

   for(int i = 1; i <= dim2; ++i)
      integrals[0] /= i;
}


// In the future, basis indices will be provided by the polygon object
void IntegralsCubeMonomial(int *dims, int deg, double *integrals)
{
   int i, d;
   int dim = dims[0];
   int m = BasisSize(deg, dim);
   INT_8 *basisIndices = (INT_8 *)malloc(dim*m*sizeof(int));
   BasisIndices(deg, dim, basisIndices);

   for(i = 0; i < m; ++i) {
      double val = 1.0;
      for(d = 0; d < dim; ++d)
         val = val/(basisIndices[i*dim+d] + 1);
      integrals[i] = val;
   }

   free(basisIndices);
}
