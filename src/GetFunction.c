/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */

#include "GetFunction.h"

#include "Quadrature.h"
#include "Phi.h"
#include "GENERAL_QUADRATURE.h"

#include <stdlib.h>
#include <stdint.h>
#include <math.h>

/* GetFunction
 * Computes an evaluates function f = b-X*W * for Newton's method.
 */
void GetFunction(const _DomainFuncs functions, const int_fast8_t *basis, const quadrature quad, double *f)
{
   const int n_nodes = quad.k;
   const int m = quad.params->num_funs;
   const int dim = quad.params->dim;
   double *xtemp = NULL;
   double *b = (double *)malloc(m*sizeof(double));
   double *phi = (double *)malloc(m*sizeof(double));

   functions.basisIntegrals(quad.params, b);
   for(int i = 0; i < m; ++i)
      f[i] = -1.0 * b[i];

   for(int j = 0; j < n_nodes; ++j)
   {
      xtemp = &quad.x[dim*j];
      functions.evalBasis(basis, xtemp, quad.params, phi);
      for(int i = 0; i < m; ++i)
         f[i] = f[i] + quad.w[j]*phi[i];
   }

   free(b);
   free(phi);
}
