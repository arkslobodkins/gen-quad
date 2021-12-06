/* Arkadijs Slobodkins
 * SMU Mathematics
 * November 2021 // Added OpenMP
 */

#include "GetJacobian.h"
#include "Basis.h"
#include <stdio.h>

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
#ifdef _OPENMP
void GetJacobianOmp(quadrature *q, CMatrix JACOBIAN)
{
   double time_start = get_cur_time();

   int dim  = q->dim;
   int deg  = q->deg;
   int rows = q->basis->numFuncs;
   int cols = q->num_nodes;
   const double *w = q->w;
   const double *x = q->x;

   int num_dims = q->num_dims;
   int dims[num_dims];
   memcpy(dims, q->dims, num_dims*sizeof(int));

   int omp_condition = OMP_CONDITION(deg, dim);
   Basis *basis = q->basis;

   #pragma omp parallel default(shared) if(omp_condition) num_threads(omp_get_max_threads())
   {
      Vector basisFuncsLoc = basis->basisOmpData.basis_funcs_omp[omp_get_thread_num()];
      Vector basisDerLoc   = basis->basisOmpData.basis_der_omp[omp_get_thread_num()];
      #pragma omp for schedule(static)
      for(int j = 0; j < cols; ++j) {
         double curNode[dim];
         for(int i = 0; i < dim; ++i)
            curNode[i] = x[dim*j+i];

         BasisFuncs(basis, curNode, basisFuncsLoc);
         BasisDer(basis, curNode, basisDerLoc);
         for(int i = 0; i < rows; ++i) {
            JACOBIAN.cid[j][i] = basisFuncsLoc.id[i];
            for(int d = 0; d < dim; ++d)
               JACOBIAN.cid[j*dim+cols+d][i] = basisDerLoc.id[d*rows+i] * w[j];
         }
      }
   } // end omp parallel

   JACOBIAN_TIME += get_cur_time() - time_start;
}
#endif

void GetJacobian(quadrature *q, CMatrix JACOBIAN)
{
   double time_start = get_cur_time();

   int dim  = q->dim;
   int rows = q->basis->numFuncs;
   int cols = q->num_nodes;
   const double *w = q->w;
   const double *x = q->x;

   Basis *basis       = q->basis;
   Vector basis_funcs = basis->basis_funcs;
   Vector basis_der   = basis->basis_der;
   for(int j = 0; j < cols; ++j)
   {
      const double *curNode = &x[dim*j];
      BasisFuncs(basis, curNode, basis_funcs);
      BasisDer(basis, curNode, basis_der);

      for(int i = 0; i < rows; ++i) {
         JACOBIAN.cid[j][i] = basis_funcs.id[i];
         for(int d = 0; d < dim; ++d)
            JACOBIAN.cid[j*dim+cols+d][i] = basis_der.id[d*rows+i] * w[j];
      }
   }

   JACOBIAN_TIME += get_cur_time() - time_start;
}

