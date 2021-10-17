/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */

#include "GaussTensor.h"
#include "Quadrature.h"

#include <assert.h>

/* NodesTensor
 * Computes and returns 2-d tensor product of nodes
 */
void NodesTensor(const_quadrature *quad1, const_quadrature *quad2, quadrature *quad_new)
{
   assert(quad1->k * quad2->k == quad_new->k);
   int n1 = quad1->k;
   int n2 = quad2->k;

   for(int i = 0; i < n1; ++i)
   {
      for(int j = 0; j < n2; ++j)
      {
         quad_new->x[2*(n2*i+j)] = quad1->x[i];
         quad_new->x[2*(n2*i+j)+1] = quad2->x[j];
      }
   }
}


/* WeightsTensor
 * Computes and returns 2-d tensor product of weights
 */
void WeightsTensor(const_quadrature *quad1, const_quadrature *quad2, quadrature *quad_new)
{
   assert(quad1->k * quad2->k == quad_new->k);
   int n1 = quad1->k;
   int n2 = quad2->k;

   for(int i = 0; i < n1; ++i)
      for(int j = 0; j < n2; ++j)
         quad_new->w[i*n2+j] = quad1->w[i] * quad2->w[j];
}

