/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */

#include "SetParams.h"

#include "BasisIndices.h"
#include "GENERAL_QUADRATURE.h"

void SetParams(int dim, int num_dims, int *dims, int deg, quadParams *params)
{
   params->dim = dim;
   params->num_dims = num_dims;
   for(int i = 0; i < num_dims; ++i) params->dims[i] = dims[i];

   params->deg = deg;
   params->num_funs = BasisSize(params->dim, deg);
}

