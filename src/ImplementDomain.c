/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */

#include "ImplementDomain.h"

#include "ComputeDomain.h"
#include "GENERAL_QUADRATURE.h"

#include <stdio.h>

void ImplementDomain(int deg, int *dims, DOMAIN_TYPE D)
{

#ifdef QUAD_DEBUG_ON
   printf("DEBUG MODE ON\n\n");
#else
   printf("DEBUG MODE OFF\n\n");
#endif

   switch(D)
   {
      case INTERVAL:
         ComputeInterval(deg);
         break;
      case CUBE:
         ComputeCube(deg, dims[0]);
         break;
      case SIMPLEX:
         ComputeSimplex(deg, dims[0]);
         break;
      case CUBESIMPLEX:
         ComputeCubeSimplex(deg, dims[0], dims[1]);
         break;
      case SIMPLEXSIMPLEX:
         ComputeSimplexSimplex(deg, dims[0], dims[1]);
         break;
      case CUBESIMPLEXSIMPLEX:
         ComputeCubeSimplexSimplex(deg, dims[0], dims[1], dims[2]);
         break;
   }
}

