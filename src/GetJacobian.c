/* Arkadijs Slobodkins
 * SMU Mathematics
 * November 2021 // Added OpenMP
 */

#include "GetJacobian.h"

#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include "get_time.h"
double JACOBIAN_TIME = 0.0;

// Computes an returns Jacobian of size num_funcs x (dim+1)*num_nodes
// of a vector X * w, where each row of X corresponds to
// a distinct polynomial of orthogonal basis, and each
// column of X and entries in w correspond to ith node and weight
// of the quadrature.
void GetJacobian(const INT_8 *basisIndices, quadrature *q, CMatrix JACOBIAN)
{
   double time_start = get_cur_time();

   int dim  = q->dim;
   int deg  = q->deg;
   int rows = q->num_funcs;
   int cols = q->num_nodes;
   const double *w = q->w;
   const double *x = q->x;

   int num_dims = q->num_dims;
   int dims[num_dims];
   memcpy(dims, q->dims, num_dims*size_int);

   #ifdef _OPENMP
   int omp_condition = OMP_CONDITION(deg, dim);
   #endif

   // achieves speedup when dim >= 3 and quadrature degree is sufficiently high
   // if dim == 3, degree should be  > 10 for speedup
   // The higher the dimension, the lower becomes the degree requirement.
   // In the future number of threads will be selected as a function of problem size.
   #ifdef _OPENMP
   #pragma omp parallel default(shared) if(omp_condition) num_threads(omp_get_max_threads())
   #endif
   {
      double *basisLoc      = (double *)malloc(SIZE_DOUBLE(rows));
      double *basisPrimeLoc = (double *)malloc(SIZE_DOUBLE(rows*dim));
      INT_8 *basisIndCopy   = (INT_8 *)malloc(dim*q->num_funcs*sizeof(INT_8));
      memcpy(basisIndCopy, basisIndices, dim*q->num_funcs*sizeof(INT_8));


      #ifdef _OPENMP
      #pragma omp for schedule(static)
      #endif
      for(int j = 0; j < cols; ++j)
      {
         double curNode[dim];
         for(int i = 0; i < dim; ++i)
            curNode[i] = x[dim*j+i];

         q->evalBasis(dims, deg, basisIndCopy, curNode, basisLoc);
         q->evalBasisDer(dims, deg, basisIndCopy, curNode, basisPrimeLoc);

         for(int i = 0; i < rows; ++i)
         {
            JACOBIAN.cid[j][i] = basisLoc[i];

            int col_index = j*dim+cols;
            for(int d = 0; d < dim; ++d)
               JACOBIAN.cid[col_index+d][i] = basisPrimeLoc[d*rows+i] * w[j];
         }
      }
      free(basisLoc);
      free(basisPrimeLoc);
      free(basisIndCopy);
   } // end omp parallel

   JACOBIAN_TIME += get_cur_time() - time_start;
}

