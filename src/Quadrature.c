/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */

#include "Quadrature.h"

#include "BasisIndices.h"
#include "BasisFunctions.h"
#include "BasisIntegrals.h"
#include "Constraints.h"
#include "Print.h"

#include <stdlib.h>
#include <assert.h>
#include <string.h>


static inline int GetNumDims(DOMAIN_TYPE D);
static inline bool QuadPosWeightsElem(const quadrature *q, int elem);

static void SetIntervalFuncs(quadrature *q);
static void SetCubeFuncs(quadrature *q);
static void SetSimplexFuncs(quadrature *q);
static void SetCubeSimplexFuncs(quadrature *q);
static void SetSimplexSimplexFuncs(quadrature *q);
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

      default:
         return NULL;
   }

   quadrature *q = (quadrature *)malloc(size_quadrature);
   memset(q, 0, sizeof(quadrature));

   q->num_nodes = n;
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
   }

   q->D = D;
   q->evalBasisMonomial = &BasisMonomial;

   q->constr           = NULL;
   q->evalBasis        = NULL;
   q->evalBasisDer     = NULL;
   q->basisIntegrals   = NULL;
   q->constr_init      = NULL;
   q->constr_realloc   = NULL;
   q->get_constr       = NULL;
   q->constr_free      = NULL;
   q->free_ptr = &quadrature_free_basic;

   q->isFullyInitialized = GQ_FALSE;

   return q;
}


quadrature *quadrature_init_full(int n, int dim, int *dims, int deg, DOMAIN_TYPE D)
{
   quadrature *q = quadrature_init_basic(n, dim, dims, deg, D);

   q->setFuncs(q);
   q->constr = q->constr_init(q->dims);
   q->get_constr(q->constr);
   q->isFullyInitialized = GQ_TRUE;
   q->free_ptr = &quadrature_free_full;

   return q;
}


// Routine that reallocates quadrature q and reinitializes its parameters,
// except for parameters corresponding to the number of dimensions
// (i.e. size of array dims remains unchanged) and DOMAIN_TYPE of q.
void quadrature_realloc(int n, int dim, int *dims, int deg, quadrature *q)
{
   q->num_nodes = n;

   Vector_realloc(n*(dim+1), &q->z);
   q->w = &q->z.id[0];
   q->x = &q->z.id[n];

   q->dim = dim;
   q->deg = deg;
   q->num_funcs = BasisSize(deg, dim);
   for(int i = 0; i < q->num_dims; ++i)
      q->dims[i] = dims[i];

   if(q->isFullyInitialized == GQ_TRUE) {
      q->constr_realloc(q->constr, dims);
      q->get_constr(q->constr);
   }

}


// Routine that reallocates q to a quadrature with n nodes.
// First n nodes are copied to the new quadrature.
// Before the routine is called, q should contain at least n nodes or more.
void quadrature_reinit(int n, quadrature *q)
{
   assert(q->num_nodes >= n);
   q->num_nodes = n;

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


quadrature *quadrature_make_full_copy(const quadrature *q)
{
   if(q->isFullyInitialized == GQ_FALSE) {
      PRINT_ERR("supplied quadrature object was not fully initialized", __LINE__, __FILE__);
      return NULL;
   }

   int k         = q->num_nodes;
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


quadrature *quadrature_without_element(quadrature *q, int i)
{
   int count, j, d;
   int num_nodes = q->num_nodes;
   int dim = q->dim;
   quadrature *q_without = quadrature_init_full(num_nodes-1, q->dim, q->dims, q->deg, q->D);

   for(count = 0, j = 0; j < num_nodes; ++j) {
      if(j == i) continue;
      q_without->w[count] = q->w[j];
      for(d = 0; d < dim; ++d)
         q_without->x[count*dim+d] = q->x[j*dim+d];
      ++count;
   }

   return q_without;
}


// Assignment operator
void quadrature_assign(const quadrature *q1, quadrature *q2)
{
   assert(q2->num_nodes == q1->num_nodes);
   int k   = q1->num_nodes;
   int dim = q1->dim;

   // assign only nodes and weights, other fields remain unchanged
   memcpy( q2->w, q1->w, SIZE_DOUBLE(k) );
   memcpy( q2->x, q1->x, SIZE_DOUBLE(k*dim) );
}


void quadrature_remove_element(int index, quadrature *q)
{
   assert(index <= q->num_nodes-1);

   int k     = q->num_nodes;
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
   assert( (q.num_nodes * (q.dim+1)) == v.len );

   int k = q.num_nodes;
   int dim = q.dim;

   memcpy( &v.id[0], q.w, SIZE_DOUBLE(k) );
   memcpy( &v.id[k], q.x, SIZE_DOUBLE(k*dim) );
}


void quadrature_get_elem(const quadrature *q, int i, Vector v)
{
   assert(v.len = q->dim);
   int dim = q->dim;
   memcpy( &v.id[0], &q->w[i], SIZE_DOUBLE(1) );
   memcpy( &v.id[1], &q->x[i*dim], SIZE_DOUBLE(dim) );
}


void vector_to_quadrature(const Vector v, quadrature q)
{
   assert( (q.num_nodes * (q.dim+1)) == v.len );

   int k = q.num_nodes;
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
   if(q->isFullyInitialized == GQ_TRUE)
   {
      PRINT_ERR("full destructor should be called", __LINE__, __FILE__);
      return;
   }
   if(q == NULL) return;

   Vector_free(q->z);
   if(q->dims != NULL) { free(q->dims); q->dims = NULL; }
   free(q); q = NULL;
}


static void quadrature_free_full(quadrature *q)
{
   if(q->isFullyInitialized == GQ_FALSE)
   {
      PRINT_ERR("basic destructor should be called", __LINE__, __FILE__);
      return;
   }
   if(q == NULL) return;

   Vector_free(q->z);
   if(q->dims != NULL) { free(q->dims); q->dims = NULL; }
   if(q->isFullyInitialized) q->constr_free(q->constr);
   if(q->dims != NULL) { free(q->dims); q->dims = NULL; }
   free(q); q = NULL;
}


bool QuadOnTheBoundary(const quadrature *q, int elem)
{
   if(q->isFullyInitialized == GQ_FALSE) {
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


bool QuadInDomain(const quadrature *q)
{
   if(q->isFullyInitialized == GQ_FALSE) {
      PRINT_ERR("constraints have not been initialized", __LINE__, __FILE__);
      return false;
   }

   int dim = q->dim;
   RMatrix M = q->constr->M;
   Vector b = q->constr->b;

   for(int i = 0; i < q->num_nodes; ++i)
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


bool QuadInDomainElem(const quadrature *q, int elem)
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


bool QuadPosWeights(const quadrature *q)
{
   const quadrature _q = *q;
   for(int i = 0; i < _q.num_nodes; ++i)
      if(_q.w[i] < 0)
         return false;

   return true;
}


static inline bool QuadPosWeightsElem(const quadrature *q, int elem)
{
   return is_greater_than_zero(q->w[elem]);
}


bool QuadInConstraint(const quadrature *q)
{
   if(q->constr == NULL) {
      PRINT_ERR("constraints have not been initialized", __LINE__, __FILE__);
      return false;
   }

   return QuadInDomain(q) & QuadPosWeights(q);
}


bool QuadInConstraintElem(const quadrature *q, int elem)
{
   if(q->constr == NULL) {
      PRINT_ERR("constraints have not been initialized", __LINE__, __FILE__);
      return false;
   }

   return QuadInDomainElem(q, elem) & QuadPosWeightsElem(q, elem);
}


double QuadTestIntegral(const quadrature *q)
{
   if(q->constr == NULL) {
      PRINT_ERR("constraints have not been initialized", __LINE__, __FILE__);
      return false;
   }

   int k = q->num_nodes;
   int dim = q->dim;
   int *dims = q->dims;
   int deg = q->deg;
   double *x = q->x;
   double *w = q->w;

   int num_funcs = q->num_funcs;
   INT_8 *basis_id = (INT_8 *)malloc(num_funcs*dim *sizeof(INT_8));
   BasisIndices(deg, dim, basis_id);

   Vector orth_basis = Vector_init(num_funcs);
   Vector IQuad      = Vector_init(num_funcs);
   Vector IExact     = Vector_init(num_funcs);
   Vector res_arr    = Vector_init(num_funcs);

   q->basisIntegrals(dims, deg, IExact.id);

   // approximate integrals of basis functions
   for(int i = 0; i < k; ++i)
   {
      q->evalBasis(dims, deg, basis_id, &x[dim*i], orth_basis.id);
      for(int j = 0; j < num_funcs; ++j)
         IQuad.id[j] += orth_basis.id[j] * w[i];
   }

   for(int j = 0; j < num_funcs; ++j) res_arr.id[j] = fabs(IQuad.id[j]-IExact.id[j]);
   double res = V_ScaledTwoNorm(res_arr);

   Vector_free(IExact);
   Vector_free(IQuad);
   Vector_free(orth_basis);
   Vector_free(res_arr);
   free(basis_id);
   return res;
}


double QuadTestIntegralMonomial(const quadrature *q)
{
   if(q->constr == NULL) {
      PRINT_ERR("constraints have not been initialized", __LINE__, __FILE__);
      return false;
   }

   int k = q->num_nodes;
   int dim = q->dim;
   int *dims = q->dims;
   int deg = q->deg;
   double *x = q->x;
   double *w = q->w;

   int num_funcs = q->num_funcs;
   INT_8 *basis_id = (INT_8 *)malloc(num_funcs*dim *sizeof(INT_8));
   BasisIndices(deg, dim, basis_id);

   Vector orth_basis = Vector_init(num_funcs);
   Vector IQuad      = Vector_init(num_funcs);
   Vector IExact     = Vector_init(num_funcs);
   Vector res_arr    = Vector_init(num_funcs);

   q->basisIntegralsMonomial(dims, deg, IExact.id);

   // approximate integrals of basis functions
   for(int i = 0; i < k; ++i)
   {
      q->evalBasisMonomial(dim, deg, basis_id, &x[dim*i], orth_basis.id);
      for(int j = 0; j < num_funcs; ++j)
         IQuad.id[j] += orth_basis.id[j] * w[i];
   }

   for(int j = 0; j < num_funcs; ++j) res_arr.id[j] = fabs(IQuad.id[j]-IExact.id[j]);
   double res = V_ScaledTwoNorm(res_arr);

   Vector_free(IExact);
   Vector_free(IQuad);
   Vector_free(orth_basis);
   Vector_free(res_arr);
   free(basis_id);
   return res;
}


bool QuadEqnOnTheBoundary(const quadrature *q, int elem, int eqn)
{
   if(q->constr == NULL || q->isFullyInitialized == GQ_FALSE) {
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
   q->evalBasis        = &CubeBasisFuncs;
   q->evalBasisDer     = &CubeBasisFuncsDer;
   q->basisIntegrals   = &BasisIntegralsCube;

   q->constr_init      = &constraints_interval_init;
   q->constr_realloc   = &constraints_interval_realloc;
   q->get_constr       = &get_constraints_interval;
   q->constr_free      = &constraints_interval_free;
}


static void SetCubeFuncs(quadrature *q)
{
   q->evalBasis        = &CubeBasisFuncs;
   q->evalBasisDer     = &CubeBasisFuncsDer;
   q->basisIntegrals   = &BasisIntegralsCube;
   q->basisIntegralsMonomial = &IntegralsCubeMonomial; // Will add to other polygons as well

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


static inline int GetNumDims(DOMAIN_TYPE D)
{
   switch(D)
   {
      case INTERVAL:           return ONE;
      case CUBE:               return ONE;
      case SIMPLEX:            return ONE;
      case CUBESIMPLEX:        return TWO;
      case SIMPLEXSIMPLEX:     return TWO;
      default:                 return -1;
   }
}


