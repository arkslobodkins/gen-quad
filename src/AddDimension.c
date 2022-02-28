
#include "AddDimension.h"
#include "GeneralGaussTensor.h"

#include <assert.h>

void AddLineFirst(const quadrature *q1D, const quadrature *quad_prev, quadrature *quad_new)
{
   assert(q1D->dim == 1);
   assert(q1D->dim + quad_prev->dim == quad_new->dim);
   assert(q1D->num_nodes * quad_prev->num_nodes == quad_new->num_nodes);

   const int n1D        = q1D->num_nodes;
   const double *x1D    = q1D->x;
   const int n_prev     = quad_prev->num_nodes;
   const double *x_prev = quad_prev->x;
   const int dim        = quad_new->dim;
   const int sizePrev   = dim*n_prev;
   double *x_new        = quad_new->x;

   // compute tensor product for nodes
   for(int i = 0; i < n1D; ++i)
      for(int j = 0; j < n_prev; ++j)
         x_new[i*sizePrev+j*dim] = x1D[i];

   for(int i = 0; i < n1D; ++i)
      for(int j = 0; j < n_prev; ++j)
         for(int d = 1; d < dim; ++d)
            x_new[i*sizePrev+j*dim+d] = x_prev[j*(dim-1)+d-1];

   WeightsTensor2D(q1D, quad_prev, quad_new);
}


void AddLineSimplex(const quadrature *q1D, const quadrature *quad_prev, quadrature *quad_new)
{
   assert(q1D->dim == 1);
   assert(q1D->dim + quad_prev->dim == quad_new->dim);
   assert(q1D->num_nodes * quad_prev->num_nodes == quad_new->num_nodes);

   const int dim = quad_new->dim;
   const int n_new = quad_new->num_nodes;
   double *x_new = quad_new->x;
   double *w_new = quad_new->w;

   AddLineFirst(q1D, quad_prev, quad_new);
   WeightsTensor2D(q1D, quad_prev, quad_new);

   // apply one level of Duffy Transformation
   for(int i = 0; i < n_new; ++i)
      for(int d = 1; d < dim; ++d)
         w_new[i] *= x_new[i*dim];

   for(int i = 0; i < n_new; ++i)
      for(int d = 1; d < dim; ++d)
         x_new[i*dim+d] *= x_new[i*dim];

}


void GeneralDuffy(quadrature *q)
{
   const int dim = q->dim;
   const int N = q->num_nodes;
   double *x = q->x;
   double *w = q->w;

   for(int i = 0; i < N; ++i)
      for(int d = 1; d < dim; ++d)
      {
         x[i*dim+d] *= x[i*dim+d-1];
         w[i] *= x[i*dim+d-1];
      }
}


void MixedTensor(const quadrature *q1, const quadrature *q2, quadrature *q_tp)
{
   assert(q1->dim + q2->dim == q_tp->dim);
   assert(q1->num_nodes * q2->num_nodes == q_tp->num_nodes);

   const int n1     = q1->num_nodes;
   const int dim1   = q1->dim;
   const double *x1 = q1->x;
   const int n2     = q2->num_nodes;
   const int dim2   = q2->dim;
   const double *x2 = q2->x;
   const int dim_tp = q_tp->dim;
   double *xtp      = q_tp->x;

   for(int i = 0; i < n1; ++i)
   {
      for(int j = 0; j < n2; ++j)
      {
         int id = i*n2*dim_tp + j*dim_tp;
         for(int d1 = 0; d1 < dim1; ++d1)
            xtp[id+d1] = x1[i*dim1+d1];

         int count = 0;
         for(int d2 = dim1; d2 < dim_tp; ++d2)
         {
            xtp[id+d2] = x2[j*dim2+count];
            ++count;
         }
      }
   }
   WeightsTensor2D(q1, q2, q_tp);
}
