/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */

#include "BasisIntegrals.h"

#include "AddDimension.h"
#include "BasisIndices.h"
#include "Gauss_Lib/Jacobi.h"
#include "Phi.h"
#include "Quadrature.h"
#include "GENERAL_QUADRATURE.h"

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

/* Routines that compute analytical integrals of all orthogonal basis functions
 * for respective domain and parameters.
 */
void IntegralsCube(const quadParams *params, double *integrals)
{
   const int m = params->num_funs;
   memset(integrals, 0, m*sizeof(integrals[0]));
   integrals[0] = 1.0;
}


void IntegralsSimplex(const quadParams *params, double *integrals)
{
   const int m = params->num_funs;
   const int dim1 = params->dims[0];

   memset(integrals, 0, m*sizeof(integrals[0]));
   integrals[0] = 1;
   for(int i = 1; i <= dim1; ++i)
      integrals[0] /= i;

}


void IntegralsCubeSimplex(const quadParams *params, double *integrals)
{
   const int m = params->num_funs;
   const int dim2 = params->dims[1];

   memset(integrals, 0, m*sizeof(integrals[0]));
   integrals[0] = 1;
   for(int i = 1; i <= dim2; ++i)
      integrals[0] /= i;
}


void IntegralsSimplexSimplex(const quadParams *params, double *integrals)
{
   const int m = params->num_funs;
   const int dim1 = params->dims[0];
   const int dim2 = params->dims[1];

   memset(integrals, 0, m*sizeof(integrals[0]));
   integrals[0] = 1;
   for(int i = 1; i <= dim1; ++i)
      integrals[0] /= i;

   for(int i = 1; i <= dim2; ++i)
      integrals[0] /= i;

}


void IntegralsCubeSimplexSimplex(const quadParams *params, double *integrals)
{
   const int m = params->num_funs;
   const int dim2 = params->dims[0];
   const int dim3 = params->dims[1];

   memset(integrals, 0, m*sizeof(integrals[0]));
   integrals[0] = 1;
   for(int i = 1; i <= dim2; ++i)
      integrals[0] /= i;

   for(int i = 1; i <= dim3; ++i)
      integrals[0] /= i;

}

