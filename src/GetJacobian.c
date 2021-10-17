/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */

#include "GetJacobian.h"

#include "GENERAL_QUADRATURE.h"
#include <stdlib.h>
#include <stdint.h>

/* GetJacobian
 * Computes an returns Jacobian of size m x 3*k
 * of a vector X * w, where each row of X corresponds to
 * a distinct polynomial of orthogonal basis, and each
 * column of X and entries in w correspond to ith node and weight
 * of the quadrature
 */
void GetJacobian(const _DomainFuncs functions, const int_fast8_t *basis, const quadrature quad, double *jacobian)
{
   int dim = quad.params->dim;
   int n_nodes = quad.k;
   int n_rows = quad.params->num_funs;
   int n_cols = n_nodes*(dim+1);
   double *curNode = NULL;
   double *phi = (double *)malloc(n_rows*sizeof(double));
   double *phiPrime = (double *)malloc(n_rows*dim*sizeof(double));

   // compute analytical derivatives with respect to all weights and node components
   for(int j = 0; j < n_nodes; ++j)
   {
      curNode = &quad.x[dim*j];
      functions.evalBasisDer(basis, &curNode[0], quad.params, &phiPrime[0]);
      functions.evalBasis(basis, &curNode[0], quad.params, &phi[0]);

      for(int i = 0; i < n_rows; ++i)
      {
         jacobian[ ij2(i, j, n_cols) ] = phi[i];

         int row_id;
         int col_id;
         int index = n_nodes+j*dim;
         for(int d = 0; d < dim; ++d)
         {
            row_id = i;
            col_id = index+d;
            jacobian[ ij2(row_id, col_id, n_cols) ] = phiPrime[i + d*n_rows] * quad.w[j];
         }
      }
   }

   free(phi);
   free(phiPrime);
}

