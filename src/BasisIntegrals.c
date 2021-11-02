/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */

#include "BasisIntegrals.h"

#include "BasisIndices.h"
#include "Quadrature.h"
#include "GENERAL_QUADRATURE.h"

#include <stdlib.h>
#include <string.h>

// Routines that compute analytical integrals of all orthogonal
// basis functions for respective domains and parameters.

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


void BasisIntegralsCubeSimplexSimplex(int *dims, int deg, double *integrals)
{
   int dim = dims[0]+dims[1]+dims[1];
   int dim2 = dims[1];
   int dim3 = dims[1];
   int m = BasisSize(deg, dim);

   memset(integrals, 0, SIZE_DOUBLE(m));
   integrals[0] = 1;
   for(int i = 1; i <= dim2; ++i)
      integrals[0] /= i;

   for(int i = 1; i <= dim3; ++i)
      integrals[0] /= i;
}

