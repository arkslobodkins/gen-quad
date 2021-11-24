/* Arkadijs Slobodkins
 * SMU Mathematics
 * November 2021 // Added OpenMP
 */

#include "GetJacobian.h"
#include "GENERAL_QUADRATURE.h"

#include <stdlib.h>
#include <string.h>
#include <omp.h>


// Computes an returns Jacobian of size num_funcs x (dim+1)*num_nodes
// of a vector X * w, where each row of X corresponds to
// a distinct polynomial of orthogonal basis, and each
// column of X and entries in w correspond to ith node and weight
// of the quadrature.
void GetJacobian(const INT_8 *basisIndices, quadrature *q, CMatrix JACOBIAN)
{
   int dim  = q->dim;
   int deg  = q->deg;
   int rows = q->num_funcs;
   int cols = q->num_nodes;
   const double *w = q->w;
   const double *x = q->x;

   int num_dims = q->num_dims;
   int dims[num_dims];
   memcpy(dims, q->dims, num_dims*size_int);

   bool_enum OMP_CONDITION;
   if(dim == 3  && deg >= 12)      OMP_CONDITION = GQ_TRUE;
   else if(dim == 4  && deg  >= 6) OMP_CONDITION = GQ_TRUE;
   else if(dim == 5  && deg  >= 5) OMP_CONDITION = GQ_TRUE;
   else if(dim > 5)                OMP_CONDITION = GQ_TRUE;
   else                            OMP_CONDITION = GQ_FALSE;

   // achieves speedup when dim >= 3 and quadrature degree is sufficiently high
   // if dim == 3, degree should be  > 10 for speedup
   // The higher the dimension, the lower becomes the degree requirement.
   #pragma omp parallel default(shared) if(OMP_CONDITION) num_threads(omp_get_max_threads())
   {
      double *basisLoc      = (double *)malloc(SIZE_DOUBLE(rows));
      double *basisPrimeLoc = (double *)malloc(SIZE_DOUBLE(rows*dim));
      INT_8 *basisIndCopy   = (INT_8 *)malloc(dim*q->num_funcs*sizeof(INT_8));
      memcpy(basisIndCopy, basisIndices, dim*q->num_funcs*sizeof(INT_8));

      #pragma omp for schedule(static)
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
}

