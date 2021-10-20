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

static inline int GetNumDims(DOMAIN_TYPE D);

static void SetIntervalFuncs(quadrature *q);
static void SetCubeFuncs(quadrature *q);
static void SetSimplexFuncs(quadrature *q);
static void SetCubeSimplexFuncs(quadrature *q);
static void SetSimplexSimplexFuncs(quadrature *q);
static void SetCubeSimplexSimplexFuncs(quadrature *q);

static bool QuadInDomain(const_quadrature *q);
static bool QuadInDomainElem(const_quadrature *q, int elem);

static bool QuadPosWeights(const_quadrature *q);
static inline bool QuadPosWeightsElem(const_quadrature *q, int elem);

static bool QuadInConstraint(const_quadrature *q);
static bool QuadInConstraintElem(const_quadrature *quad, int elem);

static bool QuadOnTheBoundary(const_quadrature *q, int elem);
static bool QuadEqnOnTheBoundary(const_quadrature *q, int elem, int eqn);

static void QuadSetFull_Constr(quadrature *q);
static double QuadTestIntegral(const_quadrature *q);


quadrature *quadrature_init(int n, int dim, int *dims, int deg, DOMAIN_TYPE D)
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

   int num_dims = GetNumDims(D);
   q->dims      = (int *)malloc(num_dims*size_int);
   for(int i = 0; i < num_dims; ++i)
      q->dims[i] = dims[i];
   q->dim       = dim;
   q->num_dims  = num_dims;
   q->deg       = deg;
   q->num_funcs = BasisSize(dim, deg);

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
   q->eqnOnTheBoundary = NULL;
   q->onTheBoundary    = NULL;
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
   q->onTheBoundary    = &QuadOnTheBoundary;
   q->eqnOnTheBoundary = &QuadEqnOnTheBoundary;

   q->setFuncs(q);
   q->cons = q->constr_init(q->dims);
   q->get_constr(q->cons);
   q->setFuncsFlag = 1;

   QuadSetFull_Constr(q);
}


// Routine that reallocates quadrature q and reinitializes its parameters,
// except for parameters corresponding to the number of dimensions
// (i.e. size of array dims remains unchanged) and DOMAIN_TYPE of q.
void quadrature_realloc(int n, int dim, int *dims, int deg, quadrature *q)
{
   q->k = n;

   q->z = (double *)realloc( q->z, SIZE_DOUBLE(n*(dim+1)) );
   q->w = &q->z[0];
   q->x = &q->z[n];

   q->dim = dim;
   q->deg = deg;
   q->num_funcs = BasisSize(dim, deg);
   for(int i = 0; i < q->num_dims; ++i)
      q->dims[i] = dims[i];
   q->setFuncs(q);
}


// Routine that reallocates q to a quadrature with n nodes.
// First n nodes are copied to the new quadrature.
// Before the routine is called, q should contain at least n nodes or more.
void quadrature_reinit(int n, quadrature *q)
{
   assert(q->k >= n);
   q->k = n;

   int dim = q->dim;
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
   int k   = q1->k;
   int dim = q1->dim;

   // assign only nodes and weights, other fields remain unchanged
   memcpy( q2->w, q1->w, SIZE_DOUBLE(k) );
   memcpy( q2->x, q1->x, SIZE_DOUBLE(k*dim) );
}


void quadrature_remove_element(int index, quadrature *q)
{
   assert(index <= q->k-1);

   int k     = q->k;
   int dim   = q->dim;
   double *w = q->w;
   double *x = q->x;

   for(int i = index; i < k-1; ++i)
      w[i] = w[i+1];

   for(int i = index; i < k-1; ++i)
      for(int d = 0; d < dim; ++d)
         x[i*dim+d] = x[(i+1)*dim+d];

   quadrature_reinit(k-1, q);
}


void quadrature_to_vector(const quadrature q, Vector v)
{
   assert( (q.k * (q.dim+1)) == v.len );

   int k = q.k;
   int dim = q.dim;

   memcpy( &v.id[0], q.w, SIZE_DOUBLE(k) );
   memcpy( &v.id[k], q.x, SIZE_DOUBLE(k*dim) );
}


void quadrature_get_elem(const_quadrature *q, int i, Vector v)
{
   int dim = q->dim;
   memcpy( &v.id[0], &q->w[i], SIZE_DOUBLE(1) );
   memcpy( &v.id[1], &q->x[i*dim], SIZE_DOUBLE(dim) );
}


void vector_to_quadrature(const Vector v, quadrature q)
{
   assert( (q.k * (q.dim+1)) == v.len );

   int k = q.k;
   int dim = q.dim;

   memcpy( q.w, v.id, SIZE_DOUBLE(k) );
   memcpy( q.x, &v.id[k], SIZE_DOUBLE(k*dim) );
}


void quadrature_free(quadrature *q)
{
   if(q == NULL) return;

   if(q->z != NULL) {
      free(q->z); q->z = NULL; }

   if(q->dims != NULL) {
      free(q->dims); q->dims = NULL; }


   if(q->setFuncsFlag)
   {
      q->constr_free(q->cons);
      Matrix_free(q->FULL_A);
      Vector_free(q->FULL_b);
   }
   if(q->dims != NULL) {
      free(q->dims); q->dims = NULL;
   }


   free(q); q = NULL;
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


static bool QuadOnTheBoundary(const_quadrature *q, int elem)
{
   if(q->cons == NULL || q->setFuncsFlag == 0)
   {
      PRINT_ERR("constraints have not been initialized", __LINE__, __FILE__);
      return false;
   }
   int dim = q->dim;
   int node_index = elem*dim;
   Vector b = q->cons->b;
   Matrix A = q->cons->M;
   int rows = A.rows;
   int cols = A.cols;
   double tol = pow(10, -12);

   for(int i = 0; i < rows; ++i)
   {
      double b_elem = 0.0;
      for(int d = 0; d < cols; ++d)
         b_elem += A.id[i][d] * q->x[node_index+d];

      if( fabs(b_elem - b.id[i]) <= tol )
         return true;
   }

   return false;
}


static bool QuadEqnOnTheBoundary(const_quadrature *q, int elem, int eqn)
{
   if(q->cons == NULL || q->setFuncsFlag == 0)
   {
      PRINT_ERR("constraints have not been initialized", __LINE__, __FILE__);
      return false;
   }
   int dim = q->dim;
   int node_index = elem*dim;
   Vector b = q->cons->b;
   Matrix A = q->cons->M;
   int cols = A.cols;
   double tol = pow(10, -12);

   double b_elem = 0.0;
   for(int d = 0; d < cols; ++d)
      b_elem += A.id[eqn][d] * q->x[node_index+d];

   if( fabs(b_elem - b.id[eqn]) <= tol )
      return true;

   return false;

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


static void SetIntervalFuncs(quadrature *q)
{
   q->evalBasis        = &PhiCube;
   q->evalBasisDer     = &PhiPrimeCube;
   q->basisIntegrals   = &IntegralsCube;
   q->constr_init      = &constraints_interval_init;
   q->constr_realloc   = &constraints_interval_realloc;
   q->get_constr       = &get_constraints_interval;
   q->constr_free      = &constraints_interval_free;
}


static void SetCubeFuncs(quadrature *q)
{
   q->evalBasis        = &PhiCube;
   q->evalBasisDer     = &PhiPrimeCube;
   q->basisIntegrals   = &IntegralsCube;
   q->constr_init      = &constraints_cube_init;
   q->constr_realloc   = &constraints_cube_realloc;
   q->get_constr       = &get_constraints_cube;
   q->constr_free      = &constraints_cube_free;
}


static void SetSimplexFuncs(quadrature *q)
{
   q->evalBasis        = &PhiSimplex;
   q->evalBasisDer     = &PhiPrimeSimplex;
   q->basisIntegrals   = &IntegralsSimplex;
   q->constr_init      = &constraints_simplex_init;
   q->constr_realloc   = &constraints_simplex_realloc;
   q->get_constr       = &get_constraints_simplex;
   q->constr_free      = &constraints_simplex_free;
}


static void SetCubeSimplexFuncs(quadrature *q)
{
   q->evalBasis        = &PhiCubeSimplex;
   q->evalBasisDer     = &PhiPrimeCubeSimplex;
   q->basisIntegrals   = &IntegralsCubeSimplex;
   q->constr_init      = &constraints_cubesimplex_init;
   q->constr_realloc   = &constraints_cubesimplex_realloc;
   q->get_constr       = &get_constraints_cubesimplex;
   q->constr_free      = &constraints_cubesimplex_free;
}


static void SetSimplexSimplexFuncs(quadrature *q)
{
   q->evalBasis        = &PhiSimplexSimplex;
   q->evalBasisDer     = &PhiPrimeSimplexSimplex;
   q->basisIntegrals   = &IntegralsSimplexSimplex;
   q->constr_init      = &constraints_simplexsimplex_init;
   q->constr_realloc   = &constraints_simplexsimplex_realloc;
   q->get_constr       = &get_constraints_simplexsimplex;
   q->constr_free      = &constraints_simplexsimplex_free;
}


static void SetCubeSimplexSimplexFuncs(quadrature *q)
{
   q->evalBasis        = &PhiCubeSimplexSimplex;
   q->evalBasisDer     = &PhiPrimeCubeSimplexSimplex;
   q->basisIntegrals   = &IntegralsCubeSimplexSimplex;
   q->constr_init      = &constraints_cubesimplexsimplex_init;
   q->constr_realloc   = &constraints_cubesimplexsimplex_realloc;
   q->get_constr       = &get_constraints_cubesimplexsimplex;
   q->constr_free      = &constraints_cubesimplexsimplex_free;
}


static bool QuadInDomain(const_quadrature *q)
{
   if(q->cons == NULL)
   {
      PRINT_ERR("constraints have not been initialized", __LINE__, __FILE__);
      return false;
   }

   int dim = q->dim;
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


static bool QuadInDomainElem(const_quadrature *q, int elem)
{
   if(q->cons == NULL)
   {
      PRINT_ERR("constraints have not been initialized", __LINE__, __FILE__);
      return false;
   }

   int dim = q->dim;
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


static bool QuadPosWeights(const_quadrature *q)
{
   const_quadrature _q = *q;
   for(int i = 0; i < _q.k; ++i)
      if(_q.w[i] < 0)
         return false;

   return true;
}


static inline bool QuadPosWeightsElem(const_quadrature *q, int elem)
{
   return q->w[elem] > 0 ? true : false;
}


static bool QuadInConstraint(const_quadrature *q)
{
   return QuadInDomain(q) & QuadPosWeights(q);
}


static bool QuadInConstraintElem(const_quadrature *quad, int elem)
{
   return QuadInDomainElem(quad, elem) & QuadPosWeightsElem(quad, elem);
}


static double QuadTestIntegral(const_quadrature *q)
{
   int num_funcs   = q->num_funcs;
   int dim = q->dim;
   int deg = q->deg;
   int_fast8_t *basis_id = (int_fast8_t *)malloc(num_funcs*dim *sizeof(int_fast8_t));
   BasisIndices(deg, dim, basis_id);

   Vector res_loc = Vector_init(num_funcs);
   Vector phi = Vector_init(num_funcs);
   Vector In = Vector_init(num_funcs);
   Vector Ie = Vector_init(num_funcs);

   q->basisIntegrals((int *)q->dims, q->deg, Ie.id);

   // approximate integrals of basis functions
   for(int i = 0; i < q->k; ++i)
   {
      q->evalBasis((int *)q->dims, q->deg, basis_id, &q->x[dim*i], phi.id);
      for(int j = 0; j < num_funcs; ++j)
         In.id[j] += phi.id[j] * q->w[i];
   }

   for(int j = 0; j < num_funcs; ++j) res_loc.id[j] = fabs(In.id[j]-Ie.id[j]);
   double res = V_ScaledTwoNorm(res_loc);

   Vector_free(Ie);
   Vector_free(In);
   Vector_free(phi);
   Vector_free(res_loc);
   free(basis_id);
   return res;
}
