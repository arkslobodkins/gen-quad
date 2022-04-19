/* Arkadijs Slobodkins
 * SMU Mathematics
 * Copyright 2022
 */


#include "GetFunction.h"
#include "Quadrature.h"
#include "Basis.h"
#include "Print.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <omp.h>

#include "get_time.h"
double FUNCTION_TIME = 0.0;

// Computes an evaluates function f = F(X)*W - b for Newton's method.
#ifdef _OPENMP
void GetFunctionOmp(quadrature *q, Vector f)
{
   assert(f.len == q->basis->numFuncs);
   assert(q->omp_threads == f.ompLen);
   double time_start = get_cur_time();

   int nthreads = f.ompLen;
   int dim = q->dim;
   int deg = q->deg;
   int len = q->basis->numFuncs;
   int nodes = q->num_nodes;
   const double *w = q->w;
   Tensor2D x = DoubleToTensor2D(nodes, dim, q->x);

   Basis **basis = q->basisOmp;
   BasisIntegrals(basis[0], basis[0]->integrals);
   Vector_Assign(basis[0]->integrals, f);
   VScale(-1.0, f);

   int omp_condition = OMP_CONDITION(deg, dim);
   #pragma omp parallel if(omp_condition) default(shared) num_threads(nthreads)
   {
      double *fLoc = f.ompId[omp_get_thread_num()];
      memset(fLoc, 0, SIZE_DOUBLE(len));
      Basis *basisLoc = basis[omp_get_thread_num()];
      Vector basisFuncsLoc = basisLoc->functions;

      #pragma omp for schedule(static)
      for(int i = 0; i < nodes; ++i)
      {
         double curNode[dim];
         for(int d = 0; d < dim; ++d)
            curNode[d] = TID2(x, i, d);

         BasisFuncs(basisLoc, curNode, basisFuncsLoc);
         double_daxpy(len, w[i], basisFuncsLoc.id, fLoc);
      }
      #pragma omp critical(UpdateFunction)
      double_daxpy(len, 1.0, fLoc, f.id);
   } // end omp parallel
   FUNCTION_TIME += get_cur_time() - time_start;
}
#endif


// Computes an evaluates function f = F(X)*W - b for Newton's method.
void GetFunction(quadrature *q, Vector f)
{
   assert(f.len == q->basis->numFuncs);
   double time_start = get_cur_time();

   int dim = q->dim;
   int nodes = q->num_nodes;
   const double *w = q->w;
   const double *x = q->x;
   Vector functions = q->basis->functions;
   Vector integrals = q->basis->integrals;

   BasisIntegrals(q->basis, integrals);
   Vector_Assign(integrals, f);
   VScale(-1.0, f);
   for(int i = 0; i < nodes; ++i)
   {
      const double *curNode = &x[dim*i];
      BasisFuncs(q->basis, curNode, functions);
      Vector_daxpy(w[i], functions, f);
   }
   FUNCTION_TIME += get_cur_time() - time_start;
}


void TestResidual(quadrature *q, const char* str)
{
   Vector f = Vector_init(q->basis->numFuncs);
   GetFunction(q, f);
   double norm = V_InfNorm(f);
   PrintDouble(norm, str);
   Vector_free(f);
}
