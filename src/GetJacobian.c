/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */

#include "GetJacobian.h"

#include "GENERAL_QUADRATURE.h"
#include <stdlib.h>
#include <stdint.h>

// Computes an returns Jacobian of size m x (dim+1)*k
// of a vector X * w, where each row of X corresponds to
// a distinct polynomial of orthogonal basis, and each
// column of X and entries in w correspond to ith node and weight
// of the quadrature
void GetJacobian(const int_fast8_t *basis, quadrature *q, CMatrix JACOBIAN)
{
   int dim          = q->dim;
   int n_nodes      = q->k;
   int n_rows       = q->num_funcs;
   double *phi      = (double *)malloc(n_rows*sizeof(double));
   double *phiPrime = (double *)malloc(n_rows*dim*sizeof(double));

   // compute analytical derivatives with respect to all weights and node components
   for(int j = 0; j < n_nodes; ++j)
   {
      double *curNode = &q->x[dim*j];
      q->evalBasisDer(q->dims, q->deg, basis, &curNode[0], &phiPrime[0]);
      q->evalBasis(q->dims, q->deg, basis, &curNode[0], &phi[0]);

      for(int i = 0; i < n_rows; ++i)
      {
         JACOBIAN.cid[j][i] = phi[i];

         int index = n_nodes+j*dim;
         for(int d = 0; d < dim; ++d)
            JACOBIAN.cid[index+d][i] = phiPrime[i + d*n_rows] * q->w[j];
      }
   }

   free(phi);
   free(phiPrime);
}

