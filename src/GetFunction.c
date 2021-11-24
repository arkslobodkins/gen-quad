/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */

#include "GetFunction.h"
#include "Quadrature.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

// Computes an evaluates function f = F(X)*W - b for Newton's method.
void GetFunction(const INT_8 *basisIndices, quadrature *q, Vector f)
{
   assert(f.len = q->num_funcs);

   int dim   = q->dim;
   int deg   = q->deg;
   int len   = q->num_funcs;
   int nodes = q->num_nodes;
   const double *w = q->w;
   const double *x = q->x;

   int num_dims = q->num_dims;
   int dims[num_dims];
   memcpy(dims, q->dims, num_dims*size_int);

   double *b = (double *)malloc(len*sizeof(double));
   q->basisIntegrals(q->dims, deg, b);
   for(int i = 0; i < len; ++i)
      f.id[i] = -1.0 * b[i];
   free(b);

   bool_enum OMP_CONDITION;
   if(dim == 3  && deg >= 12)      OMP_CONDITION = GQ_TRUE;
   else if(dim == 4  && deg  >= 6) OMP_CONDITION = GQ_TRUE;
   else if(dim == 5  && deg  >= 5) OMP_CONDITION = GQ_TRUE;
   else if(dim > 5)                OMP_CONDITION = GQ_TRUE;
   else                            OMP_CONDITION = GQ_FALSE;

   // achieves speedup when dim >= 3 and quadrature degree is sufficiently high
   // if dim == 3, degree should be  > 10 for speedup
   // The higher the dimension, the lower becomes the degree requirement.
   #pragma omp parallel if(OMP_CONDITION) default(shared) num_threads(omp_get_max_threads())
   {
      double *basisLoc   = (double *)malloc(SIZE_DOUBLE(len));
      double *fLoc     = (double *)malloc(SIZE_DOUBLE(len));
      memset(fLoc, 0, SIZE_DOUBLE(len));
      INT_8 *basisIndCopy = (INT_8 *)malloc(dim*q->num_funcs*sizeof(INT_8));
      memcpy(basisIndCopy, basisIndices, dim*q->num_funcs*sizeof(INT_8));

      #pragma omp for schedule(static)
      for(int j = 0; j < nodes; ++j)
      {
         double curNode[dim];
         for(int i = 0; i < dim; ++i)
            curNode[i] = x[dim*j+i];

         q->evalBasis(dims, deg, basisIndCopy, curNode, basisLoc);
         for(int i = 0; i < len; ++i)
            fLoc[i] += basisLoc[i] * w[j];
      }
      #pragma omp critical(UpdateFunction)
      for(int i = 0; i < len; ++i)
         f.id[i] += fLoc[i];

      free(basisIndCopy);
      free(fLoc);
      free(basisLoc);
   } // end omp parallel

}
