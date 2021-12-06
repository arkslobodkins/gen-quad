/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */

#include "GetFunction.h"
#include "Quadrature.h"
#include "Basis.h"

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
   assert(f.len = q->basis->numFuncs);
   double time_start = get_cur_time();

   int dim   = q->dim;
   int deg   = q->deg;
   int len   = q->basis->numFuncs;
   int nodes = q->num_nodes;
   const double *w = q->w;
   const double *x = q->x;

   int num_dims = q->num_dims;
   int dims[num_dims];
   memcpy(dims, q->dims, num_dims*sizeof(int));

   Basis **basis = q->basisOmp;
   BasisIntegrals(basis[0], basis[0]->integrals);
   for(int i = 0; i < len; ++i)
      f.id[i] = -1.0 * basis[0]->integrals.id[i];

   int omp_condition = OMP_CONDITION(deg, dim);
   #pragma omp parallel if(omp_condition) default(shared) num_threads(omp_get_max_threads())
   {
      double *fLoc = f.ompId[omp_get_thread_num()];
      memset(fLoc, 0, SIZE_DOUBLE(len));
      Basis *basisLoc = basis[omp_get_thread_num()];
      Vector basisFuncsLoc = basisLoc->functions;

      #pragma omp for schedule(static)
      for(int j = 0; j < nodes; ++j)
      {
         double curNode[dim];
         for(int i = 0; i < dim; ++i)
            curNode[i] = x[dim*j+i];

         BasisFuncs(basisLoc, curNode, basisFuncsLoc);
         for(int i = 0; i < len; ++i)
            fLoc[i] += basisFuncsLoc.id[i] * w[j];
      }
      #pragma omp critical(UpdateFunction)
      for(int i = 0; i < len; ++i)
         f.id[i] += fLoc[i];
   } // end omp parallel

   FUNCTION_TIME += get_cur_time() - time_start;
}
#endif

// Computes an evaluates function f = F(X)*W - b for Newton's method.
void GetFunction(quadrature *q, Vector f)
{
   assert(f.len = q->basis->numFuncs);
   double time_start = get_cur_time();

   int dim   = q->dim;
   int nodes = q->num_nodes;
   const double *w = q->w;
   const double *x = q->x;
   Vector functions = q->basis->functions;
   Vector integrals = q->basis->integrals;

   int num_dims = q->num_dims;
   int dims[num_dims];
   memcpy(dims, q->dims, num_dims*sizeof(int));

   BasisIntegrals(q->basis, integrals);
   for(int i = 0; i < integrals.len; ++i)
      f.id[i] = -1.0 * integrals.id[i];

   for(int j = 0; j < nodes; ++j)
   {
      const double *curNode = &x[dim*j];
      BasisFuncs(q->basis, curNode, functions);
      for(int i = 0; i < functions.len; ++i)
         f.id[i] += functions.id[i] * w[j];
   }
   FUNCTION_TIME += get_cur_time() - time_start;
}
