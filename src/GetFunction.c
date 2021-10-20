/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */

#include "GetFunction.h"

#include "Quadrature.h"
#include "BasisFunctions.h"
#include "GENERAL_QUADRATURE.h"

#include <stdlib.h>
#include <stdint.h>
#include <math.h>

// Computes an evaluates function f = b-X*W * for Newton's method.
void GetFunction(const int_fast8_t *basis, quadrature *q, double *f)
{
   int n_nodes   = q->k;
   int m         = q->num_funcs;
   int dim       = q->dim;
   double *b     = (double *)malloc(m*sizeof(double));
   double *phi   = (double *)malloc(m*sizeof(double));

   q->basisIntegrals(q->dims, q->deg, b);
   for(int i = 0; i < m; ++i)
      f[i] = -1.0 * b[i];

   for(int j = 0; j < n_nodes; ++j)
   {
      double *xtemp = &q->x[dim*j];
      q->evalBasis(q->dims, q->deg, basis, xtemp, phi);
      for(int i = 0; i < m; ++i)
         f[i] = f[i] + q->w[j]*phi[i];
   }

   free(b);
   free(phi);
}
