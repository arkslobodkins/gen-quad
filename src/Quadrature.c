/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */

#include "Quadrature.h"
#include "SetQuadFuncs.h"
#include "SetParams.h"
#include "GENERAL_QUADRATURE.h"
#include "Print.h"

#include <stdlib.h>
#include <assert.h>
#include <string.h>


static inline int get_num_dims(DOMAIN_TYPE D)
{
   switch(D)
   {
      case INTERVAL:           return ONE;
      case CUBE:               return ONE;
      case SIMPLEX:            return ONE;
      case CUBESIMPLEX:        return TWO;
      case SIMPLEXSIMPLEX:     return TWO;
      case CUBESIMPLEXSIMPLEX: return THREE;
      default: return -1;
   }
}


quadrature *quadrature_init(int n, int dim, int *dims, int p, DOMAIN_TYPE D)
{

   // ensure DOMAIN is initialized with valid parameters
   switch(D)
   {
      case INTERVAL:
         assert(dim == 1);
         assert(dim == dims[0]);
         break;

      case CUBE:
         assert(dim == dims[0]);
         break;

      case SIMPLEX:
         assert(dim == dims[0]);
         break;

      case CUBESIMPLEX:
         assert(dim == dims[0]+dims[1]);
         break;

      case SIMPLEXSIMPLEX:
         assert(dim == dims[0]+dims[1]);
         break;

      case CUBESIMPLEXSIMPLEX:
         assert(dim == dims[0]+dims[1]+dims[2]);
         break;

      default:
         return NULL;
   }

   quadrature *q = (quadrature *)malloc(size_quadrature);
   memset(q, 0, sizeof(quadrature));

   q->k = n;
   q->z = (double *)malloc( SIZE_DOUBLE(n*(dim+1)) );
   memset( q->z, 0, SIZE_DOUBLE(n*(dim+1)) );
   q->w = &q->z[0];
   q->x = &q->z[n];

   q->setParams = SetParams;
   q->params = (quadParams *)malloc(size_quadParams);
   int num_dims = get_num_dims(D);
   q->params->dims = (int *)malloc(SIZE_INT(num_dims));
   q->setParams(dim, num_dims, dims, p, q->params);


   switch(D)
   {
      case INTERVAL:
         q->setDomain = SetIntervalFuncs;
         break;
      case CUBE:
         q->setDomain = SetCubeFuncs;
         break;
      case SIMPLEX:
         q->setDomain = SetSimplexFuncs;
         break;
      case CUBESIMPLEX:
         q->setDomain = SetCubeSimplexFuncs;
         break;
      case SIMPLEXSIMPLEX:
         q->setDomain = SetSimplexSimplexFuncs;
         break;
      case CUBESIMPLEXSIMPLEX:
         q->setDomain = SetCubeSimplexSimplexFuncs;
         break;
   }
   q->D = D;

   q->cons = NULL;
   q->domFuncs = NULL;

   q->evalBasis      = NULL;
   q->evalBasisDer   = NULL;
   q->basisIntegrals = NULL;
   q->inDomain       = NULL;
   q->constr_init    = NULL;
   q->constr_realloc = NULL;
   q->get_constr     = NULL;
   q->constr_free    = NULL;

   return q;
}


void quadSetFuncsAndConstr(quadrature *q)
{

   if(q == NULL) return;

   if(q->domFuncs != NULL)
      PRINT_ERR("Domain functions are already set", 0, __LINE__, __FILE__);

   q->domFuncs = (DomainFuncs)malloc(sizeof(_DomainFuncs));
   q->setDomain(q);

   q->cons = q->constr_init(q->params->dims);
   q->get_constr(q->cons);
}


// Routine that reallocates quadrature q and reinitializes its parameters,
// except for parameters corresponding to the number of dimensions
// (i.e. size of array dims remains unchanged) and DOMAIN_TYPE of q.
void quadrature_realloc(int n, int dim, int *dims, int p, quadrature *q)
{
   q->k = n;

   q->z = (double *)realloc( q->z, SIZE_DOUBLE(n*(dim+1)) );
   q->w = &q->z[0];
   q->x = &q->z[n];

   q->setParams(dim, q->params->num_dims, dims, p, q->params);
   q->setDomain(q);
}


// Routine that reallocates q to a quadrature with n nodes.
// First n nodes are copied to the new quadrature.
// Before the routine is called, q should contain at least n nodes or more.
void quadrature_reinit(int n, quadrature *q)
{
   assert(q->k >= n);
   q->k = n;

   int dim = q->params->dim;
   double *x = (double *)malloc(SIZE_DOUBLE(n*dim));
   double *w = (double *)malloc(SIZE_DOUBLE(n));
   // temporarily store nodes and weights
   memcpy( w, q->w, SIZE_DOUBLE(n) );
   memcpy( x, q->x, SIZE_DOUBLE(n*dim) );

   q->z = (double *)realloc( q->z, SIZE_DOUBLE(n * (dim+1)) );
   q->w = &q->z[0];
   q->x = &q->z[n];

   // copy back
   memcpy( q->w, w, SIZE_DOUBLE(n) );
   memcpy( q->x, x, SIZE_DOUBLE(n*dim) );

   free(x);
   free(w);
}


// Assignment operator
void quadrature_assign(const quadrature q1, quadrature q2)
{
   assert(q2.k == q1.k);
   int k = q1.k;
   int dim = q1.params->dim;

   // assign only nodes and weights, other fields remain unchanged
   memcpy( q2.w, q1.w, SIZE_DOUBLE(k) );
   memcpy( q2.x, q1.x, SIZE_DOUBLE(k*dim) );
}


void quadrature_remove_element(int index, quadrature *q)
{
   assert(index <= q->k-1);

   int k = q->k;
   int dim = q->params->dim;

   for(int i = index; i < k-1; ++i)
      q->w[i] = q->w[i+1];

   for(int i = index; i < k-1; ++i)
      for(int d = 0; d < dim; ++d)
         q->x[i*dim+d] = q->x[(i+1)*dim+d];

   quadrature_reinit(k-1, q);
}


void quadrature_to_vector(const quadrature q, Vector v)
{
   assert( (q.k * (q.params->dim+1)) == v.len );

   int k = q.k;
   int dim = q.params->dim;

   memcpy( &v.id[0], q.w, SIZE_DOUBLE(k) );
   memcpy( &v.id[k], q.x, SIZE_DOUBLE(k*dim) );
}


void vector_to_quadrature(const Vector v, quadrature q)
{
   assert( (q.k * (q.params->dim+1)) == v.len );

   int k = q.k;
   int dim = q.params->dim;

   memcpy( q.w, v.id, SIZE_DOUBLE(k) );
   memcpy( q.x, &v.id[k], SIZE_DOUBLE(k*dim) );
}


void quadrature_free(quadrature *q)
{
   if(q == NULL) return;

   if(q->z != NULL) {
      free(q->z); q->z = NULL; }

   if(q->params->dims != NULL) {
      free(q->params->dims); q->params->dims = NULL; }

   if(q->params != NULL) {
      free(q->params); q->params = NULL; }

   if(q->domFuncs != NULL) {
      free(q->domFuncs); q->domFuncs = NULL; }

   if(q->cons != NULL) {
      q->constr_free(q->cons); q->cons = NULL; }

   free(q); q = NULL;
}
