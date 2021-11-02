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
static inline bool QuadPosWeightsElem(const_quadrature *q, int elem);

static void SetIntervalFuncs(quadrature *q);
static void SetCubeFuncs(quadrature *q);
static void SetSimplexFuncs(quadrature *q);
static void SetCubeSimplexFuncs(quadrature *q);
static void SetSimplexSimplexFuncs(quadrature *q);
static void SetCubeSimplexSimplexFuncs(quadrature *q);
static void quadrature_free_basic(quadrature *q);
static void quadrature_free_full(quadrature *q);



quadrature *quadrature_init_basic(int n, int dim, int *dims, int deg, DOMAIN_TYPE D)
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
   q->z = Vector_init(n*(dim+1));
   memset( q->z.id, 0, SIZE_DOUBLE(n*(dim+1)) );
   q->w = &q->z.id[0];
   q->x = &q->z.id[n];

   int num_dims = GetNumDims(D);
   q->dims      = (int *)malloc(num_dims*size_int);
   for(int i = 0; i < num_dims; ++i)
      q->dims[i] = dims[i];
   q->num_dims  = num_dims;
   q->dim       = dim;
   q->deg       = deg;
   q->num_funcs = BasisSize(deg, dim);

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

   q->constr = NULL;
   q->evalBasis        = NULL;
   q->evalBasisDer     = NULL;
   q->basisIntegrals   = NULL;
   q->constr_init      = NULL;
   q->constr_realloc   = NULL;
   q->get_constr       = NULL;
   q->constr_free      = NULL;

   q->free_ptr = &quadrature_free_basic;

   return q;
}


quadrature *quadrature_init_full(int n, int dim, int *dims, int deg, DOMAIN_TYPE D)
{
   quadrature *q = quadrature_init_basic(n, dim, dims, deg, D);

   q->setFuncs(q);

   q->constr = q->constr_init(q->dims);
   q->get_constr(q->constr);
   q->setFuncsConstrFlag = 1;

   q->free_ptr = &quadrature_free_full;

   return q;
}


// Routine that reallocates quadrature q and reinitializes its parameters,
// except for parameters corresponding to the number of dimensions
// (i.e. size of array dims remains unchanged) and DOMAIN_TYPE of q.
void quadrature_realloc(int n, int dim, int *dims, int deg, quadrature *q)
{
   q->k = n;

   Vector_realloc(n*(dim+1), &q->z);
   q->w = &q->z.id[0];
   q->x = &q->z.id[n];

   q->dim = dim;
   q->deg = deg;
   q->num_funcs = BasisSize(deg, dim);
   for(int i = 0; i < q->num_dims; ++i)
      q->dims[i] = dims[i];

   if(q->setFuncsConstrFlag == 1) {
      q->constr_realloc(q->constr, dims);
      q->get_constr(q->constr);
   }

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

   Vector_realloc(n*(dim+1), &q->z);
   q->w = &q->z.id[0];
   q->x = &q->z.id[n];

   // copy back
   memcpy( q->w, w, SIZE_DOUBLE(n) );
   memcpy( q->x, x, SIZE_DOUBLE(n*dim) );

   free(x);
   free(w);
}


quadrature *quadrature_make_full_copy(const_quadrature *q)
{
   if(q->setFuncsConstrFlag == 0) {
      PRINT_ERR("supplied quadrature object was not fully initialized", __LINE__, __FILE__);
      return NULL;
   }

   int k         = q->k;
   int dim       = q->dim;
   int deg       = q->deg;
   DOMAIN_TYPE D = q->D;

   int dims[q->num_dims];
   for(int d = 0; d < q->num_dims; ++d)
      dims[d] = q->dims[d];

   quadrature *q_copy = quadrature_init_full(k, dim, dims, deg, D);
   quadrature_assign(q, q_copy);
   return q_copy;
}

// Assignment operator
void quadrature_assign(const_quadrature *q1, quadrature *q2)
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
   assert(v.len = q->dim);
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
   q->free_ptr(q);
}


static void quadrature_free_basic(quadrature *q)
{
   if(q->setFuncsConstrFlag == 1)
   {
      PRINT_ERR("full destructor should be called", __LINE__, __FILE__);
      return;
   }

   if(q == NULL) return;

   Vector_free(q->z);

   if(q->dims != NULL) {
      free(q->dims); q->dims = NULL; }

   free(q); q = NULL;
}


static void quadrature_free_full(quadrature *q)
{
   if(q->setFuncsConstrFlag == 0)
   {
      PRINT_ERR("basic destructor should be called", __LINE__, __FILE__);
      return;
   }
   if(q == NULL) return;

   Vector_free(q->z);

   if(q->dims != NULL) {
      free(q->dims); q->dims = NULL; }

   if(q->setFuncsConstrFlag) q->constr_free(q->constr);
   if(q->dims != NULL) {
      free(q->dims); q->dims = NULL;
   }

   free(q); q = NULL;
}


bool QuadOnTheBoundary(const_quadrature *q, int elem)
{
   if(q->setFuncsConstrFlag == 0) {
      PRINT_ERR("constraints have not been initialized", __LINE__, __FILE__);
      return false;
   }
   int dim = q->dim;
   int node_index = elem*dim;
   Vector b = q->constr->b;
   RMatrix A = q->constr->M;
   int rows = A.rows;
   int cols = A.cols;
   double tol = pow(10, -12);

   for(int i = 0; i < rows; ++i)
   {
      double b_elem = 0.0;
      for(int d = 0; d < cols; ++d)
         b_elem += A.rid[i][d] * q->x[node_index+d];

      if( fabs(b_elem - b.id[i]) <= tol )
         return true;
   }

   return false;
}


bool QuadInDomain(const_quadrature *q)
{
   if(q->setFuncsConstrFlag == 0) {
      PRINT_ERR("constraints have not been initialized", __LINE__, __FILE__);
      return false;
   }

   int dim = q->dim;
   RMatrix M = q->constr->M;
   Vector b = q->constr->b;

   for(int i = 0; i < q->k; ++i)
   {
      double lhs[M.rows]; memset(lhs, 0, M.rows*size_double);
      const double *X_ixdim = &q->x[i*dim];

      for(int r = 0; r < M.rows; ++r)
         for(int c = 0; c < M.cols; ++c)
            lhs[r] += M.rid[r][c] * X_ixdim[c];

      for(int r = 0; r < M.rows; ++r)
         if(lhs[r] > b.id[r])
            return false;
   }

   return true;
}


bool QuadInDomainElem(const_quadrature *q, int elem)
{
   if(q->constr == NULL) {
      PRINT_ERR("constraints have not been initialized", __LINE__, __FILE__);
      return false;
   }

   int dim = q->dim;
   RMatrix M = q->constr->M;
   Vector b = q->constr->b;
   const double *X_ixdim = &q->x[elem*dim];


   double lhs[M.rows]; memset(lhs, 0, M.rows*size_double);

   for(int r = 0; r < M.rows; ++r)
      for(int c = 0; c < M.cols; ++c)
         lhs[r] += M.rid[r][c] * X_ixdim[c];

   for(int r = 0; r < M.rows; ++r)
      if(lhs[r] > b.id[r])
         return false;

   return true;
}


bool QuadPosWeights(const_quadrature *q)
{
   const_quadrature _q = *q;
   for(int i = 0; i < _q.k; ++i)
      if(_q.w[i] < 0)
         return false;

   return true;
}


static inline bool QuadPosWeightsElem(const_quadrature *q, int elem)
{
   return is_greater_than_zero(q->w[elem]);
}


bool QuadInConstraint(const_quadrature *q)
{
   if(q->constr == NULL) {
      PRINT_ERR("constraints have not been initialized", __LINE__, __FILE__);
      return false;
   }

   return QuadInDomain(q) & QuadPosWeights(q);
}


bool QuadInConstraintElem(const_quadrature *q, int elem)
{
   if(q->constr == NULL) {
      PRINT_ERR("constraints have not been initialized", __LINE__, __FILE__);
      return false;
   }

   return QuadInDomainElem(q, elem) & QuadPosWeightsElem(q, elem);
}


double QuadTestIntegral(const_quadrature *q)
{
   if(q->constr == NULL) {
      PRINT_ERR("constraints have not been initialized", __LINE__, __FILE__);
      return false;
   }

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


bool QuadEqnOnTheBoundary(const_quadrature *q, int elem, int eqn)
{
   if(q->constr == NULL || q->setFuncsConstrFlag == 0)
   {
      PRINT_ERR("constraints have not been initialized", __LINE__, __FILE__);
      return false;
   }
   int dim = q->dim;
   int node_index = elem*dim;
   Vector b = q->constr->b;
   RMatrix A = q->constr->M;
   int cols = A.cols;

   double b_elem = 0.0;
   for(int d = 0; d < cols; ++d)
      b_elem += A.rid[eqn][d] * q->x[node_index+d];

   if( fabs(b_elem - b.id[eqn]) <= BOUND_TOL )
      return true;

   return false;

}


/*****************************************************************
\* Implementation of polymorphic behaviour for quadrature object \*
*****************************************************************/

static void SetIntervalFuncs(quadrature *q)
{
   q->evalBasis        = &BasisCube;
   q->evalBasisDer     = &BasisPrimeCube;
   q->basisIntegrals   = &BasisIntegralsCube;

   q->constr_init      = &constraints_interval_init;
   q->constr_realloc   = &constraints_interval_realloc;
   q->get_constr       = &get_constraints_interval;
   q->constr_free      = &constraints_interval_free;
}


static void SetCubeFuncs(quadrature *q)
{
   q->evalBasis        = &BasisCube;
   q->evalBasisDer     = &BasisPrimeCube;
   q->basisIntegrals   = &BasisIntegralsCube;

   q->constr_init      = &constraints_cube_init;
   q->constr_realloc   = &constraints_cube_realloc;
   q->get_constr       = &get_constraints_cube;
   q->constr_free      = &constraints_cube_free;
}


static void SetSimplexFuncs(quadrature *q)
{
   q->evalBasis        = &BasisSimplex;
   q->evalBasisDer     = &BasisPrimeSimplex;
   q->basisIntegrals   = &BasisIntegralsSimplex;

   q->constr_init      = &constraints_simplex_init;
   q->constr_realloc   = &constraints_simplex_realloc;
   q->get_constr       = &get_constraints_simplex;
   q->constr_free      = &constraints_simplex_free;
}


static void SetCubeSimplexFuncs(quadrature *q)
{
   q->evalBasis        = &BasisCubeSimplex;
   q->evalBasisDer     = &BasisPrimeCubeSimplex;
   q->basisIntegrals   = &BasisIntegralsCubeSimplex;

   q->constr_init      = &constraints_cubesimplex_init;
   q->constr_realloc   = &constraints_cubesimplex_realloc;
   q->get_constr       = &get_constraints_cubesimplex;
   q->constr_free      = &constraints_cubesimplex_free;
}


static void SetSimplexSimplexFuncs(quadrature *q)
{
   q->evalBasis        = &BasisSimplexSimplex;
   q->evalBasisDer     = &BasisPrimeSimplexSimplex;
   q->basisIntegrals   = &BasisIntegralsSimplexSimplex;

   q->constr_init      = &constraints_simplexsimplex_init;
   q->constr_realloc   = &constraints_simplexsimplex_realloc;
   q->get_constr       = &get_constraints_simplexsimplex;
   q->constr_free      = &constraints_simplexsimplex_free;
}


static void SetCubeSimplexSimplexFuncs(quadrature *q)
{
   q->evalBasis        = &BasisCubeSimplexSimplex;
   q->evalBasisDer     = &BasisPrimeCubeSimplexSimplex;
   q->basisIntegrals   = &BasisIntegralsCubeSimplexSimplex;

   q->constr_init      = &constraints_cubesimplexsimplex_init;
   q->constr_realloc   = &constraints_cubesimplexsimplex_realloc;
   q->get_constr       = &get_constraints_cubesimplexsimplex;
   q->constr_free      = &constraints_cubesimplexsimplex_free;
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
      default:                 return -1;
   }
}


