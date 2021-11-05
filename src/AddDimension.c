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
void AddLineFirst(const quadrature *q_gauss, const quadrature *quad_prev, quadrature *quad_new)
{
   assert(q_gauss->dim + quad_prev->dim == quad_new->dim);
   assert(q_gauss->k * quad_prev->k == quad_new->k);

   int dim = quad_new->dim;
   int n = q_gauss->k;
   int n_prev = quad_prev->k;
   const double *x_gauss = q_gauss->x;
   const double *x_prev = quad_prev->x;
   double *x_new = quad_new->x;

   // compute tensor product for nodes
   for(int i = 0; i < n; ++i)
   {
      for(int j = 0; j < n_prev; ++j)
      {
         int ind = i*dim*n_prev + j*dim;
         int jxdim_minus_1 = j*(dim-1);
         x_new[ind] = x_gauss[i];
         for(int d = 1; d < dim; ++d)
            x_new[ind+d] = x_prev[jxdim_minus_1+d-1];
      }
   }
   WeightsTensor2D(q_gauss, quad_prev, quad_new);
}


// AddLineLastSimplex
// Computes tensor product of unit (dim-1)-dimensional simplex
// and 1-dimensional interval [0, 1], and maps the product to the unit
// simplex of dimension dim using Duffy transformation.
void AddLineSimplex(const quadrature *q_gauss, const quadrature *quad_prev, quadrature *quad_new)
{
   assert(q_gauss->dim + quad_prev->dim == quad_new->dim);
   assert(q_gauss->k * quad_prev->k == quad_new->k);

   int dim = quad_new->dim;
   int n = q_gauss->k;
   int n_prev = quad_prev->k;
   int n_new = quad_new->k;
   const double *x_gauss = q_gauss->x;
   const double *x_prev = quad_prev->x;
   double *x_new = quad_new->x;
   double *w_new = quad_new->w;

   // compute tensor product for nodes
   for(int i = 0; i < n; ++i)
   {
      for(int j = 0; j < n_prev; ++j)
      {
         int ind = i*dim*n_prev + j*dim;
         int jxdim_minus_1 = j*(dim-1);
         x_new[ind] = x_gauss[i];
         for(int d = 1; d < dim; ++d)
            x_new[ind+d] = x_prev[jxdim_minus_1+d-1];
      }
   }

   WeightsTensor2D(q_gauss, quad_prev, quad_new);
   // apply one level of Duffy Transformation
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


void MixedTensor(const quadrature *q1, const quadrature *q2, quadrature *q_pr)
{
   assert(q1->dim + q2->dim == q_pr->dim);
   assert(q1->k * q2->k == q_pr->k);

   int n1 = q1->k;
   int n2 = q2->k;
   int dim1 = q1->dim;
   int dim2 = q2->dim;
   int dim_ss = q_pr->dim;

   for(int i = 0; i < n1; ++i)
   {
      for(int j = 0; j < n2; ++j)
      {
         int id = i*n2*dim_ss + j*dim_ss;
         for(int d1 = 0; d1 < dim1; ++d1)
            q_pr->x[id+d1] = q1->x[i*dim1+d1];

         int count = 0;
         for(int d2 = dim1; d2 < dim_ss; ++d2)
         {
            q_pr->x[id+d2] = q2->x[j*dim2+count];
            ++count;
         }
      }
   }
   WeightsTensor2D(q1, q2, q_pr);
}
