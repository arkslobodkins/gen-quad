/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */

#include "AddDimension.h"

#include "GeneralGaussTensor.h"
#include "GENERAL_QUADRATURE.h"

#include <assert.h>

// Computes tensor product of arbitrary quadrature over (dim-1)-dimensional
// domain Ω and 1-dimensional interval [0, 1]. Quadrature of dimension dim
// and degree p over interval Ω x [0,1] is generated. Stores interval quadrature nodes
// as the first coordinate.
void AddLineFirst(const_quadrature *q_gauss, const_quadrature *quad_prev, quadrature *quad_new)
{
   assert(q_gauss->dim + quad_prev->dim == quad_new->dim);
   assert(q_gauss->k * quad_prev->k == quad_new->k);

   int dim = quad_new->dim;
   int n = q_gauss->k;
   int n_prev = quad_prev->k;
   double *x_new = quad_new->x;
   const double *x_prev = quad_prev->x;

   // compute tensor product for nodes
   for(int i = 0; i < n; ++i)
   {
      for(int j = 0; j < n_prev; ++j)
      {
         int ind = i*dim*n_prev + j*dim;
         int jxdim_minus_1 = j*(dim-1);
         x_new[ind] = q_gauss->x[i];

         for(int d = 1; d < dim; ++d)
            x_new[ind+d] = x_prev[jxdim_minus_1+d-1];
      }
   }
   WeightsTensor2D((const_quadrature *)q_gauss, quad_prev, quad_new);
}



// AddLineLastSimplex
// Computes tensor product of unit (dim-1)-dimensional simplex
// and 1-dimensional interval [0, 1], and maps the product to the unit
// simplex of dimension dim using Duffy transformation.
void AddLineSimplex(const_quadrature *q_gauss, const_quadrature *quad_prev, quadrature *quad_new)
{
   assert(q_gauss->dim + quad_prev->dim == quad_new->dim);
   assert(q_gauss->k * quad_prev->k == quad_new->k);

   int dim = quad_new->dim;
   int n = q_gauss->k;
   int n_prev = quad_prev->k;
   int n_new = quad_new->k;
   double *x_new = quad_new->x;
   double *w_new = quad_new->w;
   const double *x_prev = quad_prev->x;

   // compute tensor product for nodes
   for(int i = 0; i < n; ++i)
   {
      for(int j = 0; j < n_prev; ++j)
      {
         int ind = i*dim*n_prev + j*dim;
         x_new[ind] = q_gauss->x[i];
         for(int d = 1; d < dim; ++d)
            x_new[ind+d] = x_prev[j*(dim-1)+d-1];

      }
   }

   WeightsTensor2D(q_gauss, quad_prev, quad_new);
   // apply Duffy Transformation
   for(int i = 0; i < n_new; ++i)
   {
      int ind = i*dim;
      for(int d = 1; d < dim; ++d)
      {
         x_new[ind+d] *= x_new[ind];
         w_new[i] *= x_new[ind];
      }
   }

}


// Maps nodes and weights from cube of dimension dim to
// simplex of dimension dim using generalized Duffy transformation.
void GeneralDuffy(quadrature *quad)
{
   int dim = quad->dim;
   int N = quad->k;
   double *x = quad->x;
   double *w = quad->w;

   for(int i = 0; i < N; ++i)
   {
      int ind = i*dim;
      for(int d = 1; d < dim; ++d)
      {
         x[ind+d] *= x[ind+d-1];
         w[i] *= x[ind+d-1];
      }
   }

}
