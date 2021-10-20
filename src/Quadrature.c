/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */

#include "Quadrature.h"

#include "BasisIndices.h"
#include "BasisFunctions.h"
#include "BasisIntegrals.h"
#include "Constraints.h"
#include "GENERAL_QUADRATURE.h"
#include "Print.h"

#include <stdlib.h>
#include <assert.h>
#include <string.h>

static inline bool QuadPosWeightsElem(const_quadrature *q, int elem);

static inline void QuadSetParams(int dim, int num_dims, int *dims, int deg, quadrature *q)
{
   q->params->dim = dim;
   q->params->num_dims = num_dims;
   q->params->deg = deg;
   q->params->num_funs = BasisSize(q->params->dim, deg);
   for(int i = 0; i < num_dims; ++i)
      q->params->dims[i] = dims[i];
}


static void QuadSetFull_Constr(quadrature *q)
{
   if(q->cons == NULL || q->setFuncsFlag == 0)
   {
      PRINT_ERR("constraints have not been initialized", __LINE__, __FILE__);
      return ;
   }

   Matrix A = q->cons->M;
   Vector b = q->cons->b;
   Matrix FULL_A = Matrix_init(A.rows+1, A.cols+1);
   Vector FULL_b = Vector_init(b.len+1);

   for(int i = 0; i < A.rows; ++i)
      for(int j = 0; j < A.cols; ++j)
         FULL_A.id[i+1][j+1] = A.id[i][j];
   FULL_A.id[0][0] = -1.0;
   for(int i = 0; i < b.len; ++i)
      FULL_b.id[i+1] = b.id[i];
   FULL_b.id[0] = 0.0;

   q->FULL_A = FULL_A;
   q->FULL_b = FULL_b;
}


static inline int GetNumDims(DOMAIN_TYPE D)
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


void SetIntervalFuncs(quadrature *q)
{
   q->evalBasis        = &PhiCube;
   q->evalBasisDer     = &PhiPrimeCube;
   q->basisIntegrals   = &IntegralsCube;
   q->constr_init      = &constraints_interval_init;
   q->constr_realloc   = &constraints_interval_realloc;
   q->get_constr       = &get_constraints_interval;
   q->constr_free      = &constraints_interval_free;
}


void SetCubeFuncs(quadrature *q)
{
   q->evalBasis        = &PhiCube;
   q->evalBasisDer     = &PhiPrimeCube;
   q->basisIntegrals   = &IntegralsCube;
   q->constr_init      = &constraints_cube_init;
   q->constr_realloc   = &constraints_cube_realloc;
   q->get_constr       = &get_constraints_cube;
   q->constr_free      = &constraints_cube_free;
}


void SetSimplexFuncs(quadrature *q)
{
   q->evalBasis        = &PhiSimplex;
   q->evalBasisDer     = &PhiPrimeSimplex;
   q->basisIntegrals   = &IntegralsSimplex;
   q->constr_init      = &constraints_simplex_init;
   q->constr_realloc   = &constraints_simplex_realloc;
   q->get_constr       = &get_constraints_simplex;
   q->constr_free      = &constraints_simplex_free;
}


void SetCubeSimplexFuncs(quadrature *q)
{
   q->evalBasis        = &PhiCubeSimplex;
   q->evalBasisDer     = &PhiPrimeCubeSimplex;
   q->basisIntegrals   = &IntegralsCubeSimplex;
   q->constr_init      = &constraints_cubesimplex_init;
   q->constr_realloc   = &constraints_cubesimplex_realloc;
   q->get_constr       = &get_constraints_cubesimplex;
   q->constr_free      = &constraints_cubesimplex_free;
}


void SetSimplexSimplexFuncs(quadrature *q)
{
   q->evalBasis        = &PhiSimplexSimplex;
   q->evalBasisDer     = &PhiPrimeSimplexSimplex;
   q->basisIntegrals   = &IntegralsSimplexSimplex;
   q->constr_init      = &constraints_simplexsimplex_init;
   q->constr_realloc   = &constraints_simplexsimplex_realloc;
   q->get_constr       = &get_constraints_simplexsimplex;
   q->constr_free      = &constraints_simplexsimplex_free;
}


void SetCubeSimplexSimplexFuncs(quadrature *q)
{
   q->evalBasis        = &PhiCubeSimplexSimplex;
   q->evalBasisDer     = &PhiPrimeCubeSimplexSimplex;
   q->basisIntegrals   = &IntegralsCubeSimplexSimplex;
   q->constr_init      = &constraints_cubesimplexsimplex_init;
   q->constr_realloc   = &constraints_cubesimplexsimplex_realloc;
   q->get_constr       = &get_constraints_cubesimplexsimplex;
   q->constr_free      = &constraints_cubesimplexsimplex_free;
}


bool QuadInDomain(const_quadrature *q)
{
   if(q->cons == NULL)
   {
      PRINT_ERR("constraints have not been initialized", __LINE__, __FILE__);
      return false;
   }

   int dim = q->params->dim;
   Matrix M = q->cons->M;
   Vector b = q->cons->b;

   for(int i = 0; i < q->k; ++i)
   {
      double lhs[M.rows]; memset(lhs, 0, M.rows*size_double);
      const double *X_ixdim = &q->x[i*dim];

      for(int r = 0; r < M.rows; ++r)
         for(int c = 0; c < M.cols; ++c)
            lhs[r] += M.id[r][c] * X_ixdim[c];

      for(int r = 0; r < M.rows; ++r)
         if(lhs[r] > b.id[r])
            return false;
   }

   return true;
}



bool QuadInDomainElem(const_quadrature *q, int elem)
{
   if(q->cons == NULL)
   {
      PRINT_ERR("constraints have not been initialized", __LINE__, __FILE__);
      return false;
   }

   int dim = q->params->dim;
   Matrix M = q->cons->M;
   Vector b = q->cons->b;
   const double *X_ixdim = &q->x[elem*dim];


   double lhs[M.rows]; memset(lhs, 0, M.rows*size_double);

   for(int r = 0; r < M.rows; ++r)
      for(int c = 0; c < M.cols; ++c)
         lhs[r] += M.id[r][c] * X_ixdim[c];

   for(int r = 0; r < M.rows; ++r)
      if(lhs[r] > b.id[r])
         return false;

   return true;
}


bool QuadPosWeights(const_quadrature *q)
{
   for(int i = 0; i < q->k; ++i)
      if(q->w[i] < 0)
         return false;

   return true;
}


static inline bool QuadPosWeightsElem(const_quadrature *q, int elem)
{
   return q->w[elem] > 0 ? true : false;
}


bool QuadInConstraint(const_quadrature *q)
{
   return q->inDomain(q) & q->posWeights(q);
}


bool QuadInConstraintElem(const_quadrature *quad, int elem)
{
   return quad->inDomainElem(quad, elem) & quad->posWeightsElem(quad, elem);
}


double QuadTestIntegral(const_quadrature *q)
{
   int m = q->params->num_funs;
   int dim = q->params->dim;
   int deg = q->params->deg;
   int_fast8_t *basis_id = (int_fast8_t *)malloc(m*dim *sizeof(int_fast8_t));
   BasisIndices(deg, dim, basis_id);

   Vector res_loc = Vector_init(m);
   Vector phi = Vector_init(m);
   Vector In = Vector_init(m);
   Vector Ie = Vector_init(m);

   q->basisIntegrals(q->params, Ie.id);

   // approximate integrals of basis functions
   for(int i = 0; i < q->k; ++i)
   {
      q->evalBasis(basis_id, &q->x[dim*i], q->params, phi.id);
      for(int j = 0; j < m; ++j)
         In.id[j] += phi.id[j] * q->w[i];
   }

   for(int j = 0; j < m; ++j) res_loc.id[j] = fabs(In.id[j]-Ie.id[j]);
   double res = V_ScaledTwoNorm(res_loc);

   Vector_free(Ie);
   Vector_free(In);
   Vector_free(phi);
   Vector_free(res_loc);
   free(basis_id);
   return res;
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

   q->setParams = QuadSetParams;
   q->params = (quadParams *)malloc(size_quadParams);
   int num_dims = GetNumDims(D);
   q->params->dims = (int *)malloc(SIZE_INT(num_dims));
   q->setParams(dim, num_dims, dims, p, q);


   switch(D)
   {
      case INTERVAL:
         q->setFuncs = SetIntervalFuncs;
         break;
      case CUBE:
         q->setFuncs = SetCubeFuncs;
         break;
      case SIMPLEX:
         q->setFuncs = SetSimplexFuncs;
         break;
      case CUBESIMPLEX:
         q->setFuncs = SetCubeSimplexFuncs;
         break;
      case SIMPLEXSIMPLEX:
         q->setFuncs = SetSimplexSimplexFuncs;
         break;
      case CUBESIMPLEXSIMPLEX:
         q->setFuncs = SetCubeSimplexSimplexFuncs;
         break;
   }
   q->D = D;

   q->cons = NULL;
   q->evalBasis        = NULL;
   q->evalBasisDer     = NULL;
   q->basisIntegrals   = NULL;
   q->inDomain         = NULL;
   q->inDomainElem     = NULL;
   q->inConstraint     = NULL;
   q->inConstraintElem = NULL;
   q->posWeights       = NULL;
   q->posWeightsElem   = NULL;
   q->testIntegral     = NULL;
   q->constr_init      = NULL;
   q->constr_realloc   = NULL;
   q->get_constr       = NULL;
   q->constr_free      = NULL;

   return q;
}


void quad_set_funcs_and_constr(quadrature *q)
{
   if(q == NULL) return;
   if(q->setFuncsFlag == 1)
      PRINT_ERR("Quadrature functions and constraints" "are already set", __LINE__, __FILE__);

   q->inDomain         = &QuadInDomain;
   q->inDomainElem     = &QuadInDomainElem;
   q->inConstraint     = &QuadInConstraint;
   q->inConstraintElem = &QuadInConstraintElem;
   q->posWeights       = &QuadPosWeights;
   q->posWeightsElem   = &QuadPosWeightsElem;
   q->testIntegral     = &QuadTestIntegral;

   q->setFuncs(q);
   q->cons = q->constr_init(q->params->dims);
   q->get_constr(q->cons);
   q->setFuncsFlag = 1;

   QuadSetFull_Constr(q);
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

   q->setParams(dim, q->params->num_dims, dims, p, q);
   q->setFuncs(q);
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
void quadrature_assign(const quadrature *q1, quadrature *q2)
{
   assert(q2->k == q1->k);
   int k = q1->k;
   int dim = q1->params->dim;

   // assign only nodes and weights, other fields remain unchanged
   memcpy( q2->w, q1->w, SIZE_DOUBLE(k) );
   memcpy( q2->x, q1->x, SIZE_DOUBLE(k*dim) );
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


void quadrature_get_elem(const_quadrature *q, int i, Vector v)
{
   int dim = q->params->dim;
   memcpy( &v.id[0], &q->w[i], SIZE_DOUBLE(1) );
   memcpy( &v.id[1], &q->x[i*dim], SIZE_DOUBLE(dim) );
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

   if(q->setFuncsFlag)
   {
      q->constr_free(q->cons);
      Matrix_free(q->FULL_A);
      Vector_free(q->FULL_b);
   }


   free(q); q = NULL;
}
