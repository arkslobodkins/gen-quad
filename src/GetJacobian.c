/* Arkadijs Slobodkins
 * SMU Mathematics
 * Copyright 2022
 */


#include "GetJacobian.h"
#include "LINALG.h"
#include "Basis.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include "get_time.h"
double JACOBIAN_TIME = 0.0;

#ifdef _OPENMP
void GetJacobianOmp(quadrature *q, CMatrix JACOBIAN)
{
   double time_start = get_cur_time();
   assert(q->basis->numFuncs == JACOBIAN.rows);
   assert(q->z.len == JACOBIAN.cols);

   int dim = q->dim;
   int deg = q->deg;
   int rows = q->basis->numFuncs;
   int cols = q->num_nodes;
   const double *w = q->w;
   const double *x = q->x;
   Basis **basis = q->basisOmp;

   int num_dims = q->num_dims;
   int dims[num_dims];
   memcpy(dims, q->dims, num_dims*sizeof(int));

   int omp_condition = OMP_CONDITION(deg, dim);
   #pragma omp parallel default(shared) if(omp_condition) num_threads(omp_get_max_threads())
   {
      Basis *basisLoc      = basis[omp_get_thread_num()];
      Vector basisFuncsLoc = basisLoc->functions;
      Vector basisDerLoc   = basisLoc->derivatives;
      #pragma omp for schedule(static)
      for(int j = 0; j < cols; ++j) {
         double curNode[dim];
         for(int i = 0; i < dim; ++i)
            curNode[i] = x[dim*j+i];

         BasisFuncs(basisLoc, curNode, basisFuncsLoc);
         BasisDer(basisLoc, curNode, basisDerLoc);
         VScale(w[j], basisDerLoc);

         CMatrix_LoadToColumn(j, JACOBIAN, basisFuncsLoc);
         int colId = j*dim+cols;
         for(int d = 0; d < dim; ++d)
            CMatrix_LoadToColumnRange(d*rows, colId+d, JACOBIAN, basisDerLoc);
      }
   } // end omp parallel

   JACOBIAN_TIME += get_cur_time() - time_start;
}


void GetFunctionAndJacobianOmp(quadrature *q, Vector f, CMatrix JACOBIAN)
{
   double time_start = get_cur_time();
   assert(q->basis->numFuncs == JACOBIAN.rows);
   assert(q->z.len == JACOBIAN.cols);

   int dim = q->dim;
   int deg = q->deg;
   int rows = q->basis->numFuncs;
   int cols = q->num_nodes;
   const double *w = q->w;
   const double *x = q->x;
   Basis **basis = q->basisOmp;

   int num_dims = q->num_dims;
   int dims[num_dims];
   memcpy(dims, q->dims, num_dims*sizeof(int));

   BasisIntegrals(basis[0], basis[0]->integrals);
   int len = q->basis->numFuncs;
   for(int i = 0; i < len; ++i)
      f.id[i] = -1.0 * basis[0]->integrals.id[i];

   int omp_condition = OMP_CONDITION(deg, dim);
   #pragma omp parallel default(shared) if(omp_condition) num_threads(omp_get_max_threads())
   {
      double *fLoc = f.ompId[omp_get_thread_num()];
      memset(fLoc, 0, SIZE_DOUBLE(len));
      Basis *basisLoc      = basis[omp_get_thread_num()];
      Vector basisFuncsLoc = basisLoc->functions;
      Vector basisDerLoc   = basisLoc->derivatives;
      #pragma omp for schedule(static)
      for(int j = 0; j < cols; ++j)
      {
         double curNode[dim];
         for(int i = 0; i < dim; ++i)
            curNode[i] = x[dim*j+i];

         BasisFuncs(basisLoc, curNode, basisFuncsLoc);
         double_daxpy(len, w[j], basisFuncsLoc.id, fLoc);
         BasisDer(basisLoc, curNode, basisDerLoc);
         VScale(w[j], basisDerLoc);

         CMatrix_LoadToColumn(j, JACOBIAN, basisFuncsLoc);
         int colId = j*dim+cols;
         for(int d = 0; d < dim; ++d)
            CMatrix_LoadToColumnRange(d*rows, colId+d, JACOBIAN, basisDerLoc);
      }
      #pragma omp critical(UpdateFunction)
      double_daxpy(len, 1.0, fLoc, f.id);
   } // end omp parallel

   JACOBIAN_TIME += get_cur_time() - time_start;
}
#endif


void GetJacobian(quadrature *q, CMatrix JACOBIAN)
{
   double time_start = get_cur_time();
   assert(q->basis->numFuncs == JACOBIAN.rows);
   assert(q->z.len == JACOBIAN.cols);

   int dim = q->dim;
   int rows = q->basis->numFuncs;
   int cols = q->num_nodes;
   const double *w = q->w;
   const double *x = q->x;

   Basis *basis       = q->basis;
   Vector functions   = basis->functions;
   Vector derivatives = basis->derivatives;
   for(int j = 0; j < cols; ++j)
   {
      const double *curNode = &x[dim*j];
      BasisFuncs(basis, curNode, functions);
      BasisDer(basis, curNode, derivatives);
      VScale(w[j], derivatives);

      CMatrix_LoadToColumn(j, JACOBIAN, functions);
      int colId = j*dim+cols;
      for(int d = 0; d < dim; ++d)
         CMatrix_LoadToColumnRange(d*rows, colId+d, JACOBIAN, derivatives);
   }
   JACOBIAN_TIME += get_cur_time() - time_start;
}


void GetFunctionAndJacobian(quadrature *q, Vector f, CMatrix JACOBIAN)
{
   double time_start = get_cur_time();
   assert(q->basis->numFuncs == JACOBIAN.rows);
   assert(q->z.len == JACOBIAN.cols);
   int dim = q->dim;
   int rows = q->basis->numFuncs;
   int cols = q->num_nodes;
   const double *w = q->w;
   const double *x = q->x;

   Basis *basis       = q->basis;
   Vector functions   = basis->functions;
   Vector derivatives = basis->derivatives;
   Vector integrals   = q->basis->integrals;

   BasisIntegrals(q->basis, integrals);
   for(int i = 0; i < integrals.len; ++i)
      f.id[i] = -1.0 * integrals.id[i];

   for(int j = 0; j < cols; ++j)
   {
      const double *curNode = &x[dim*j];
      BasisFuncs(basis, curNode, functions);
      Vector_daxpy(w[j], functions, f);
      CMatrix_LoadToColumn(j, JACOBIAN, functions);

      BasisDer(basis, curNode, derivatives);
      VScale(w[j], derivatives);
      int colId = j*dim+cols;
      for(int d = 0; d < dim; ++d)
         CMatrix_LoadToColumnRange(d*rows, colId+d, JACOBIAN, derivatives);
   }
   JACOBIAN_TIME += get_cur_time() - time_start;
}


void GetBasis(quadrature *q, CMatrix BasisMatrix)
{
   int dim = q->dim;
   int rows = q->basis->numFuncs;
   int cols = q->num_nodes;
   const double *x = q->x;

   assert(BasisMatrix.rows == rows);
   assert(BasisMatrix.cols == cols);

   Basis *basis = q->basis;
   Vector functions = basis->functions;
   for(int j = 0; j < cols; ++j)
   {
      const double *curNode = &x[dim*j];
      BasisFuncs(basis, curNode, functions);
      CMatrix_LoadToColumn(j, BasisMatrix, functions);
   }
}


Vector MinSingvJacobians(int n, quadrature **q)
{
   Vector min_singv = Vector_init(n);

   int rows[n];
   int cols[n];
   for(int i = 0; i < n; ++i)
   {
      rows[i] = q[i]->basis->numFuncs;
      cols[i] = q[i]->num_nodes * (q[i]->dim+1);
   }

   CMatrix Jacobians[n];
   for(int i = 0; i < n; ++i)
      Jacobians[i] = CMatrix_init(rows[i], cols[i]);

   for(int i = 0; i < n; ++i)
      GetJacobian(q[i], Jacobians[i]);
   for(int i = 0; i < n; ++i)
      min_singv.id[i] = MIN_SINGV_LAPACK(Jacobians[i]);

   for(int i = 0; i < n; ++i)
      CMatrix_free(Jacobians[i]);

   return min_singv;
}
