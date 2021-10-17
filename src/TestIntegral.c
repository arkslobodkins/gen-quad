/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */

#include "TestIntegral.h"

#include "BasisIndices.h"
#include "Vector.h"
#include "GENERAL_QUADRATURE.h"

#include <stdlib.h>
#include <stdint.h>
#include <math.h>

/* TestIntegral
 * Computes the 2-norm of the difference between the sum of all analytical
 * integrals of polynomial basis functions and their quadrature approximations.
 * The difference sum (residual) is returned.
 */
double TestIntegral(const_quadrature *q, const _DomainFuncs funcs)
{
   int m = q->params->num_funs;
   int dim = q->params->dim;
   int deg = q->params->deg;
   int_fast8_t *basis_id = (int_fast8_t *)malloc(m*dim *sizeof(int_fast8_t));
   BasisIndices(deg, dim, basis_id);

   Vector res_loc = Vector_init(m);
   Vector phi = Vector_init(m);
   Vector In = Vector_init(m);
   Vector Ie = Vector_init(m);

   funcs.basisIntegrals(q->params, Ie.id);

   // approximate integrals of basis functions
   for(int i = 0; i < q->k; ++i)
   {
      funcs.evalBasis(basis_id, &q->x[dim*i], q->params, phi.id);
      for(int j = 0; j < m; ++j)
         In.id[j] += phi.id[j] * q->w[i];
   }

   for(int j = 0; j < m; ++j) res_loc.id[j] = fabs(In.id[j]-Ie.id[j]);
   double res = V_ScaledTwoNorm(res_loc);

   Vector_free(Ie);
   Vector_free(In);
   Vector_free(phi);
   Vector_free(res_loc);
   free(basis_id);
   return res;
}

