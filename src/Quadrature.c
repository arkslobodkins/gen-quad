/* Arkadijs Slobodkins
 * SMU Mathematics
 * Copyright 2022
 */


#include "Quadrature.h"
#include "AddDimension.h"
#include "Basis.h"
#include "GeneralGaussTensor.h"
#include "Constraints.h"
#include "Print.h"

#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <time.h>

static inline int GetNumDims(DOMAIN_TYPE D);
static inline bool QuadPosWeightsElem(const quadrature *q, int elem);
static inline double QuadDistAlphaEqn(const quadrature *q, int elem, int eqn);

static void CopyBasis(const quadrature *q1, quadrature *q2);
static void SetIntervalFuncs(quadrature *q);
static void SetCubeFuncs(quadrature *q);
static void SetSimplexFuncs(quadrature *q);
static void SetCubeSimplexFuncs(quadrature *q);
static void SetSimplexSimplexFuncs(quadrature *q);
static void quadrature_free_basic(quadrature *q);
static void quadrature_free_full(quadrature *q);

static void SetCubeBasis(quadrature *q);
static void SetCubeConstr(quadrature *q);
static void SetSimplexBasis(quadrature *q);
static void SetSimplexConstr(quadrature *q);
static void SetCubeSimplexBasis(quadrature *q);
static void SetCubeSimplexConstr(quadrature *q);
static void SetSimplexSimplexBasis(quadrature *q);
static void SetSimplexSimplexConstr(quadrature *q);

quadrature* quadrature_init_basic(int n, int dim, int *dims, int deg, DOMAIN_TYPE D)
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
         assert(dims[0] >= 2);
         break;

      case SIMPLEX:
         assert(dim == dims[0]);
         assert(dims[0] >= 2);
         break;

      case CUBESIMPLEX:
         assert(dim == dims[0]+dims[1]);
         assert(dims[0] >= 1);
         assert(dims[1] >= 2);
         break;

      case SIMPLEXSIMPLEX:
         assert(dim == dims[0]+dims[1]);
         assert(dims[0] >= 2);
         assert(dims[1] >= 2);
         break;

      default:
         return NULL;
   }

   quadrature *q = (quadrature *)malloc(sizeof(quadrature));
   memset(q, 0, sizeof(quadrature));

   q->num_nodes = n;
   q->z = Vector_uninitialized(n*(dim+1));
   memset( q->z.id, 0, SIZE_DOUBLE(n*(dim+1)) );
   q->w = &q->z.id[0];
   q->x = &q->z.id[n];

   int num_dims = GetNumDims(D);
   q->dims = (int *)malloc(num_dims*sizeof(int));
   for(int i = 0; i < num_dims; ++i)
      q->dims[i] = dims[i];
   q->num_dims = num_dims;
   q->dim = dim;
   q->deg = deg;
   q->D = D;
   switch(D)
   {
      case INTERVAL:
         q->setBasisAndConstr = SetIntervalFuncs;
         break;
      case CUBE:
         q->setBasisAndConstr = SetCubeFuncs;
         q->setBasis          = SetCubeBasis;
         q->setConstr         = SetCubeConstr;
         break;
      case SIMPLEX:
         q->setBasisAndConstr = SetSimplexFuncs;
         q->setBasis          = SetSimplexBasis;
         q->setConstr         = SetSimplexConstr;
         break;
      case CUBESIMPLEX:
         q->setBasisAndConstr = SetCubeSimplexFuncs;
         q->setBasis          = SetCubeSimplexBasis;
         q->setConstr         = SetCubeSimplexConstr;
         break;
      case SIMPLEXSIMPLEX:
         q->setBasisAndConstr = SetSimplexSimplexFuncs;
         q->setBasis          = SetSimplexSimplexBasis;
         q->setConstr         = SetSimplexSimplexConstr;
         break;
   }
   q->constr           = NULL;
   q->basis            = NULL;
   q->basisOmp         = NULL;
   q->omp_threads      = 0;
   q->expIntegralExact = NULL;
   q->free_ptr = &quadrature_free_basic;
   q->isFullyInitialized = GQ_FALSE;

   return q;
}


quadrature* quadrature_init_full(int n, int dim, int *dims, int deg, DOMAIN_TYPE D)
{
   quadrature *q = quadrature_init_basic(n, dim, dims, deg, D);

   q->setBasisAndConstr(q);
   q->free_ptr = &quadrature_free_full;
   q->isFullyInitialized = GQ_TRUE;
   return q;
}


quadrature* quadrature_gauss_legendre(int deg)
{
   double p1 = 0.0, p2 = 0.0;
   int d_gauss[1] = {1};
   int n = ceil( (deg+1)/2.0 );
   quadrature *q_gauss = quadrature_init_basic(n, 1,  d_gauss, 2*n-1, INTERVAL);
   Jacobi(q_gauss->num_nodes, p1, p2, q_gauss->x, q_gauss->w);

   return q_gauss;
}


quadrature* quadrature_gauss_jacobi(int deg, double alpha, double beta)
{
   int d_gauss[1] = {1};
   int nn = ceil( (deg+1)/2.0 );
   quadrature *q_gauss = quadrature_init_basic(nn, 1,  d_gauss, 2*nn-1, INTERVAL);
   Jacobi(q_gauss->num_nodes, alpha, beta, q_gauss->x, q_gauss->w);

   // scale by a factor corresponding to G-L weights
   double c = exp(lgamma(1.0+alpha)) *
              exp(lgamma(1.0+beta)) /
              exp(lgamma(2.0+alpha+beta));
   for(int i = 0; i < q_gauss->num_nodes; ++i)
      q_gauss->w[i] *= c;

   return q_gauss;
}


quadrature* quadrature_full_cube_tensor(int deg, int dim)
{
   assert(dim > 1);
   assert(deg >= 1);

   quadrature *q1D = quadrature_gauss_legendre(deg);

   int nnt = POW_INT(q1D->num_nodes, dim);
   int d_init[1] = {dim};
   quadrature *qt = quadrature_init_full(nnt, dim, d_init, deg, CUBE);

   GeneralizedNodesTensor(q1D, qt);
   GeneralizedWeightsTensor(q1D, qt);

   quadrature_free(q1D);
   return qt;
}


quadrature* quadrature_full_simplex_tensor(int deg, int dim)
{
   assert(dim > 1);
   assert(deg >= 1);

   quadrature *qj = quadrature_gauss_jacobi(deg, 1, 0.0);
   quadrature *ql = quadrature_gauss_legendre(deg);

   int dims[1] = {2};
   quadrature *q_next = quadrature_init_full(qj->num_nodes * ql->num_nodes, 2, dims, deg, SIMPLEX);
   AddLineSimplex(qj, ql, q_next);

   quadrature *q_prev = quadrature_make_full_copy(q_next);
   for(int d = 3; d <= dim; ++d)
   {
      quadrature *q1D = quadrature_gauss_jacobi(deg, d-1, 0);

      quadrature_free(q_next);
      dims[0] = d;
      q_next = quadrature_init_full(q1D->num_nodes * q_prev->num_nodes, d, dims, deg, SIMPLEX);
      AddLineSimplex(q1D, q_prev, q_next);

      quadrature_free(q_prev);
      q_prev = quadrature_make_full_copy(q_next);

      quadrature_free(q1D);
   }

   quadrature_free(qj);
   quadrature_free(ql);
   quadrature_free(q_prev);

   return q_next;
}


void quadrature_reinit_basic(int n, int dim, int *dims, int deg, quadrature *q)
{
   if(q->isFullyInitialized == GQ_TRUE)
   {
      PRINT_ERR("routine only supports basic quadrature objects", __LINE__, __FILE__);
      return ;
   }
   q->num_nodes = n;

   Vector_realloc(n*(dim+1), &q->z);
   q->w = &q->z.id[0];
   q->x = &q->z.id[n];

   q->dim = dim;
   q->deg = deg;
   for(int i = 0; i < q->num_dims; ++i)
      q->dims[i] = dims[i];
}


void quadrature_realloc_array(int n, quadrature *q)
{
   if(q->num_nodes == n) return;

   q->num_nodes = n;
   int dim = q->dim;
   Vector_realloc(n*(dim+1), &q->z);
   q->w = &q->z.id[0];
   q->x = &q->z.id[n];
}


void quadrature_shrink_array(int n, quadrature *q)
{
   assert(q->num_nodes >= n);
   if(q->num_nodes == n) return;

   int dim = q->dim;
   double *x = (double *)malloc(SIZE_DOUBLE(n*dim));
   double *w = (double *)malloc(SIZE_DOUBLE(n));
   // temporarily store nodes and weights
   memcpy(w, q->w, SIZE_DOUBLE(n));
   memcpy(x, q->x, SIZE_DOUBLE(n*dim));

   Vector_realloc(n*(dim+1), &q->z);
   q->num_nodes = n;
   q->w = &q->z.id[0];
   q->x = &q->z.id[n];

   // copy back
   memcpy(q->w, w, SIZE_DOUBLE(n));
   memcpy(q->x, x, SIZE_DOUBLE(n*dim));

   free(x);
   free(w);
}


quadrature* quadrature_make_full_copy(const quadrature *q)
{
   if(q->isFullyInitialized == GQ_FALSE)
   {
      PRINT_ERR(STR_QUAD_NOT_FULL_INIT, __LINE__, __FILE__);
      return NULL;
   }
   int k         = q->num_nodes;
   int dim       = q->dim;
   int deg       = q->deg;
   DOMAIN_TYPE D = q->D;

   int dims[q->num_dims];
   #pragma omp simd
   for(int d = 0; d < q->num_dims; ++d)
      dims[d] = q->dims[d];

   quadrature *q_copy = quadrature_init_basic(k, dim, dims, deg, D);
   CopyBasis(q, q_copy);
   q->setConstr(q_copy);

   q_copy->expIntegralExact = q->expIntegralExact;
   q_copy->free_ptr = &quadrature_free_full;
   q_copy->isFullyInitialized = GQ_TRUE;

   quadrature_assign(q, q_copy);
   return q_copy;
}


quadrature* quadrature_without_element(const quadrature *q, int i)
{
   int count, j, d;
   int num_nodes = q->num_nodes;
   int dim = q->dim;
   quadrature *q_without = quadrature_init_full(num_nodes-1, q->dim, q->dims, q->deg, q->D);

   for(count = 0, j = 0; j < num_nodes; ++j)
   {
      if(j == i) continue;
      q_without->w[count] = q->w[j];
      #pragma omp simd
      for(d = 0; d < dim; ++d)
         q_without->x[count*dim+d] = q->x[j*dim+d];
      ++count;
   }
   return q_without;
}


void quadrature_assign(const quadrature *q1, quadrature *q2)
{
   assert(q2->num_nodes == q1->num_nodes);
   int k   = q1->num_nodes;
   int dim = q1->dim;

   // assign only nodes and weights, other fields remain unchanged
   memcpy(q2->w, q1->w, SIZE_DOUBLE(k));
   memcpy(q2->x, q1->x, SIZE_DOUBLE(k*dim));
}


void quadrature_assign_resize(const quadrature *q1, quadrature *q2)
{
   int resize = q1->num_nodes;
   quadrature_realloc_array(resize, q2);
   quadrature_assign(q1, q2);
}


void quadrature_remove_element(int index, quadrature *q)
{
   assert(index <= q->num_nodes-1);

   int k = q->num_nodes;
   int dim = q->dim;
   double *w = q->w;
   double *x = q->x;

   for(int i = index; i < k-1; ++i)
      w[i] = w[i+1];

   for(int i = index; i < k-1; ++i)
      for(int d = 0; d < dim; ++d)
         x[i*dim+d] = x[(i+1)*dim+d];

   quadrature_shrink_array(k-1, q);
}


void quadrature_to_vector(const quadrature q, Vector v)
{
   assert((q.num_nodes * (q.dim+1)) == v.len);
   int k = q.num_nodes;
   int dim = q.dim;
   memcpy(&v.id[0], q.w, SIZE_DOUBLE(k));
   memcpy(&v.id[k], q.x, SIZE_DOUBLE(k*dim));
}


void quadrature_get_elem(const quadrature *q, int i, Vector v)
{
   assert(v.len == q->dim+1);
   int dim = q->dim;
   memcpy(&v.id[0], &q->w[i], SIZE_DOUBLE(1));
   memcpy(&v.id[1], qnode(q, i), SIZE_DOUBLE(dim));
}


void vector_to_quadrature(const Vector v, quadrature q)
{
   assert( (q.num_nodes * (q.dim+1)) == v.len );
   int k = q.num_nodes;
   int dim = q.dim;
   memcpy(q.w, v.id, SIZE_DOUBLE(k));
   memcpy(q.x, &v.id[k], SIZE_DOUBLE(k*dim));
}


void quadrature_fill_random(quadrature *q)
{
   double randVal = (double)rand() / (double)RAND_MAX;
   srand(time(0) + 100.*randVal);
   for(int i = 0; i < q->z.len; ++i)
      q->z.id[i] = randVal;
}


void quadrature_free(quadrature *q)
{
   q->free_ptr(q);
}


static void quadrature_free_basic(quadrature *q)
{
   if(q->isFullyInitialized == GQ_TRUE)
   {
      PRINT_ERR("full destructor must be called instead", __LINE__, __FILE__);
      return;
   }
   if(q == NULL) return;

   Vector_free(q->z);
   if(q->dims != NULL) { free(q->dims); q->dims = NULL; }
   free(q);
}


static void quadrature_free_full(quadrature *q)
{
   if(q->isFullyInitialized == GQ_FALSE)
   {
      PRINT_ERR("basic destructor must be called instead", __LINE__, __FILE__);
      return;
   }
   if(q == NULL) return;

   Vector_free(q->z);
   if(q->dims != NULL) { free(q->dims); q->dims = NULL; }
   if(q->dims != NULL) { free(q->dims); q->dims = NULL; }
   if(q->constr != NULL) { constraints_free(q->constr); q->constr = NULL; }
   if(q->basis != NULL) { BasisFree(q->basis); q->basis = NULL;}
   free(q);
}


#ifdef _OPENMP
void QuadAllocBasisOmp(quadrature *q, int num_threads)
{
   if(q->isFullyInitialized == GQ_FALSE)
   {
      PRINT_ERR(STR_QUAD_NOT_FULL_INIT, __LINE__, __FILE__);
      return;
   }
   q->basisOmp = (Basis **)malloc(num_threads*sizeof(Basis *));
   for(int i = 0; i < num_threads; ++i)
      q->basisOmp[i] = BasisInit(q->basis->params, q->basis->interface);
   q->omp_threads = num_threads;
}


void QuadFreeBasisOmp(quadrature *q)
{
   for(int i = 0; i < q->omp_threads; ++i)
      BasisFree(q->basisOmp[i]);
   free(q->basisOmp); q->basisOmp = NULL;
}
#endif


////////////////////////////////////////////////////////////////////
bool QuadInDomainElem(const quadrature *q, int elem)
{
   if(q->isFullyInitialized == GQ_FALSE)
   {
      PRINT_ERR(STR_QUAD_NOT_FULL_INIT, __LINE__, __FILE__);
      return false;
   }

   RMatrix M = q->constr->M;
   Vector b = q->constr->b;
   Vector vnode = DToVec(q->dim, qnode(q, elem));
   StaticVectorInit(M.rows, lhs);

   RMatVec(M, vnode, lhs);
   for(int r = 0; r < M.rows; ++r)
      if(lhs.id[r] > b.id[r])
         return false;

   return true;
}


bool QuadInDomainElemEps(const quadrature *q, int elem)
{
   if(q->isFullyInitialized == GQ_FALSE)
   {
      PRINT_ERR(STR_QUAD_NOT_FULL_INIT, __LINE__, __FILE__);
      return false;
   }

   double eps = POW_DOUBLE(10, -12);
   RMatrix M = q->constr->M;
   Vector b = q->constr->b;
   Vector vnode = DToVec(q->dim, qnode(q, elem));
   StaticVectorInit(M.rows, lhs);

   RMatVec(M, vnode, lhs);
   for(int r = 0; r < M.rows; ++r)
      if(lhs.id[r] > b.id[r] + eps)
         return false;

   return true;
}


bool QuadInDomain(const quadrature *q)
{
   if(q->isFullyInitialized == GQ_FALSE)
   {
      PRINT_ERR(STR_QUAD_NOT_FULL_INIT, __LINE__, __FILE__);
      return false;
   }

   for(int i = 0; i < q->num_nodes; ++i)
      if(!QuadInDomainElem(q, i))
            return false;
   return true;
}


bool QuadInDomainEps(const quadrature *q)
{
   if(q->isFullyInitialized == GQ_FALSE)
   {
      PRINT_ERR(STR_QUAD_NOT_FULL_INIT, __LINE__, __FILE__);
      return false;
   }

   for(int i = 0; i < q->num_nodes; ++i)
      if(!QuadInDomainElemEps(q, i))
         return false;
   return true;
}


bool QuadInDomainEqnElemEps(const quadrature *q, int elem, int eqn)
{
   if(q->isFullyInitialized == GQ_FALSE)
   {
      PRINT_ERR(STR_QUAD_NOT_FULL_INIT, __LINE__, __FILE__);
      return false;
   }
   double eps = POW_DOUBLE(10, -12);
   RMatrix M = q->constr->M;
   Vector b = q->constr->b;
   assert(eqn >= 0 && eqn < M.rows);

   double *node = qnode(q, elem);
   double lhs   = DDot(M.cols, M.rid[eqn], node);
   if(lhs > b.id[eqn] + eps) return false;
   return true;
}


////////////////////////////////////////////////////////////////////
bool QuadPosWeights(const quadrature *q)
{
   const quadrature _q = *q;
   for(int i = 0; i < _q.num_nodes; ++i)
      if(_q.w[i] < 0)
         return false;

   return true;
}


bool QuadPosWeightsEps(const quadrature *q)
{
   const quadrature _q = *q;
   double eps = POW_DOUBLE(10, -14);
   for(int i = 0; i < _q.num_nodes; ++i)
      if(_q.w[i] + eps < 0)
         return false;

   return true;
}


////////////////////////////////////////////////////////////////////
bool QuadInConstraintElem(const quadrature *q, int elem)
{
   if(q->isFullyInitialized == GQ_FALSE)
   {
      PRINT_ERR(STR_QUAD_NOT_FULL_INIT, __LINE__, __FILE__);
      return false;
   }
   return QuadInDomainElem(q, elem) & QuadPosWeightsElem(q, elem);
}


bool QuadInConstraint(const quadrature *q)
{
   if(q->isFullyInitialized == GQ_FALSE)
   {
      PRINT_ERR(STR_QUAD_NOT_FULL_INIT, __LINE__, __FILE__);
      return false;
   }
   return QuadInDomain(q) & QuadPosWeights(q);
}


bool QuadInConstraintEps(const quadrature *q)
{
   if(q->isFullyInitialized == GQ_FALSE)
   {
      PRINT_ERR(STR_QUAD_NOT_FULL_INIT, __LINE__, __FILE__);
      return false;
   }
   return QuadInDomainEps(q) & QuadPosWeightsEps(q);
}


////////////////////////////////////////////////////////////////////
bool QuadOnTheBoundary(const quadrature *q, int elem)
{
   if(q->isFullyInitialized == GQ_FALSE)
   {
      PRINT_ERR(STR_QUAD_NOT_FULL_INIT, __LINE__, __FILE__);
      return false;
   }

   Vector b = q->constr->b;
   RMatrix A = q->constr->M;
   Vector vnode = DToVec(q->dim, qnode(q, elem));
   double tol = POW_DOUBLE(10.0, -12);

   for(int i = 0; i < A.rows; ++i)
   {
      double b_elem = RDotRow(i, A, vnode);
      if( fabs(b_elem - b.id[i]) <= tol )
         return true;
   }
   return false;
}


bool QuadEqnOnTheBoundary(const quadrature *q, int elem, int eqn)
{
   if(q->isFullyInitialized == GQ_FALSE)
   {
      PRINT_ERR(STR_QUAD_NOT_FULL_INIT, __LINE__, __FILE__);
      return false;
   }
   double *node = qnode(q, elem);
   Vector b = q->constr->b;
   RMatrix A = q->constr->M;

   double b_elem = DDot(A.cols, A.rid[eqn], node);
   if(fabs(b_elem - b.id[eqn]) <= BOUND_TOL)
      return true;

   return false;
}


double QuadDistFromTheBoundaryElem(const quadrature *q, int elem)
{
   if(q->isFullyInitialized == GQ_FALSE)
   {
      PRINT_ERR(STR_QUAD_NOT_FULL_INIT, __LINE__, __FILE__);
      return -1.0;
   }
   double *node = qnode(q, elem);
   Vector b = q->constr->b;
   RMatrix A = q->constr->M;

   double minDist = 0.0;
   for(int i = 0; i < A.rows; ++i)
   {
      double b_elem = DDot(A.cols, A.rid[i], node);
      if(i == 0) minDist = fabs( fabs(b_elem) - b.id[0] );
      else       minDist = MIN( fabs( fabs(b_elem) - b.id[i] ), minDist );
   }
   return minDist;
}


double QuadMinDistFromTheBoundary(const quadrature *q)
{
   if(q->isFullyInitialized == GQ_FALSE)
   {
      PRINT_ERR(STR_QUAD_NOT_FULL_INIT, __LINE__, __FILE__);
      return -1.0;
   }
   double minDist = QuadDistFromTheBoundaryElem(q, 0);
   for(int i = 1; i < q->num_nodes; ++i)
      minDist = MIN(minDist, QuadDistFromTheBoundaryElem(q, i));
   return minDist;
}

double QuadDistFromTheBoundaryTwoNorm(const quadrature *q)
{
   if(q->isFullyInitialized == GQ_FALSE)
   {
      PRINT_ERR(STR_QUAD_NOT_FULL_INIT, __LINE__, __FILE__);
      return -1.0;
   }

   Vector dist = Vector_init(q->num_nodes);
   for(int i = 0; i < q->num_nodes; ++i)
      dist.id[i] = QuadDistFromTheBoundaryElem(q, i);

   double TwoNorm = V_TwoNorm(dist);

   Vector_free(dist);
   return TwoNorm;
}


static inline bool QuadPosWeightsElem(const quadrature *q, int elem)
{
   return is_greater_than_zero(q->w[elem]);
}


static inline double QuadDistAlphaEqn(const quadrature *q, int elem, int eqn)
{
   assert(q->dim == q->constr->M.cols);

   int dim = q->dim;
   RMatrix C = q->constr->M;
   Vector b = q->constr->b;
   double weight = qweight(q, elem);
   double *node = qnode(q, elem);

   double rhs   = DDot(dim, C.rid[eqn], node);
   double alpha = (b.id[eqn] - rhs) / weight;
   return alpha;
}

double QuadDistAlphaElem(const quadrature *q, int elem)
{
   int numEqns = q->constr->M.rows;
   double minAlpha = QuadDistAlphaEqn(q, elem, 0);

   for(int j = 1; j < numEqns; ++j)
      minAlpha = MIN(minAlpha, QuadDistAlphaEqn(q, elem, j));
   return minAlpha;
}


VMin QuadMinAlpha(const quadrature *q)
{
   if(QuadInDomain(q) == false) PRINT_ERR(STR_INV_INPUT, __LINE__, __FILE__);
   int numNodes = q->num_nodes;
   Vector minAlpha = Vector_init(numNodes);

   for(int i = 0; i < numNodes; ++i)
      minAlpha.id[i] = QuadDistAlphaElem(q, i);

   VMin min = VectorMin(minAlpha);
   Vector_free(minAlpha);

   return min;
}


bool QuadIsEqual(const quadrature *a, const quadrature *b)
{
   if(a->dim != b->dim)             return false;
   if(a->num_dims != b->num_dims)   return false;
   if(a->num_nodes != b->num_nodes) return false;
   if(a->D != b->D)                 return false;

   for(int i = 0; i < a->num_dims; ++i)
      if(a->dims[i] != b->dims[i])  return false;

   for(int i = 0; i < a->z.len; ++i)
         if(!CompareDouble(a->z.id[i], b->z.id[i]))
            return false;
   return true;
}


double QuadTestIntegral(const quadrature *q, BASIS_TYPE btype)
{
   if(q->isFullyInitialized == GQ_FALSE)
   {
      PRINT_ERR(STR_QUAD_NOT_FULL_INIT, __LINE__, __FILE__);
      return -1.0;
   }
   int numFuncs = q->basis->numFuncs;

   Vector IQuad   = Vector_init(numFuncs);
   Vector res_arr = Vector_init(numFuncs);
   Vector integrals = q->basis->integrals;
   BasisFuncsPtr basisPtr = NULL;
   switch(btype)
   {
      case orthogonal:
         BasisIntegrals(q->basis, integrals);
         basisPtr = &BasisFuncs;
         break;
      case monomial:
         BasisIntegralsMonomial(q->basis, integrals);
         basisPtr = &BasisMonomial;
         break;
   }
   Vector functions = q->basis->functions;

   for(int i = 0; i < q->num_nodes; ++i)
   {
      basisPtr(q->basis, qnode(q, i), functions);
      for(int j = 0; j < numFuncs; ++j)
         IQuad.id[j] += functions.id[j] * qweight(q, i);
   }

   #pragma omp simd
   for(int j = 0; j < numFuncs; ++j) res_arr.id[j] = fabs(IQuad.id[j]-integrals.id[j]);
   double res = V_InfNorm(res_arr);

   Vector_free(IQuad);
   Vector_free(res_arr);
   return res;
}


double QuadTestIntegralExp(const quadrature *q)
{
   if(q->isFullyInitialized == GQ_FALSE)
   {
      PRINT_ERR(STR_QUAD_NOT_FULL_INIT, __LINE__, __FILE__);
      return -1.0;
   }

   double IExact = q->expIntegralExact(q);
   double IQuad = 0.0;
   for(int i = 0; i < q->num_nodes; ++i)
      IQuad += expNDim(q->dim, qnode(q, i)) * qweight(q, i);

   double res = fabs(IQuad - IExact);
   return res;
}



static void CopyBasis(const quadrature *q1, quadrature *q2)
{
   assert(q1->D == q2->D);
   q2->basis = BasisCopy(q1->basis);
}

static double ExpIntegralExactCube(const quadrature *q)
{
   return expIntegralNDimCube(q->dim);
}

static double ExpIntegralExactSimplex(const quadrature *q)
{
   return expIntegralNDimSimplex(q->dim);
}

static double ExpIntegralExactCubeSimplex(const quadrature *q)
{
   return expIntegralNDimCube(q->dims[0])*expIntegralNDimSimplex(q->dims[1]);
}

static double ExpIntegralExactSimplexSimplex(const quadrature *q)
{
   return expIntegralNDimSimplex(q->dims[0])*expIntegralNDimSimplex(q->dims[1]);
}



/*****************************************************************
\* Implementation of polymorphic behaviour for quadrature object \*
*****************************************************************/

static void SetIntervalFuncs(quadrature *q)
{
   dimParamsInterval params;
   constrInterface constrinterface = SetIntervalConstrInterface();
   q->constr = constraints_init((void *)&params, &constrinterface);
   constraints_get(q->constr);
}


static void SetCubeFuncs(quadrature *q)
{
   assert(q->D == CUBE);
   q->expIntegralExact = &ExpIntegralExactCube;
   SetCubeBasis(q);
   SetCubeConstr(q);
}


static void SetSimplexFuncs(quadrature *q)
{
   assert(q->D == SIMPLEX);
   q->expIntegralExact = &ExpIntegralExactSimplex;
   SetSimplexBasis(q);
   SetSimplexConstr(q);
}


static void SetCubeSimplexFuncs(quadrature *q)
{
   assert(q->D == CUBESIMPLEX);
   q->expIntegralExact = &ExpIntegralExactCubeSimplex;
   SetCubeSimplexBasis(q);
   SetCubeSimplexConstr(q);
}


static void SetSimplexSimplexFuncs(quadrature *q)
{
   assert(q->D == SIMPLEXSIMPLEX);
   q->expIntegralExact = &ExpIntegralExactSimplexSimplex;
   SetSimplexSimplexBasis(q);
   SetSimplexSimplexConstr(q);
}


static void SetCubeBasis(quadrature *q)
{
   assert(q->D == CUBE);
   CubeParams cubeParams;
   cubeParams.deg = q->deg;
   cubeParams.dim = q->dim;
   CubeParams *params = &cubeParams;
   BasisInterface interface = SetCubeBasisInterface();
   Basis *cube = BasisInit((void *)params, &interface);
   q->basis = cube;
}


static void SetCubeConstr(quadrature *q)
{
   assert(q->D == CUBE);
   dimParamsCube dimCube = {q->dim};
   constrInterface constrinterface = SetCubeConstrInterface();
   q->constr = constraints_init((void *)&dimCube, &constrinterface);
   constraints_get(q->constr);
}


static void SetSimplexBasis(quadrature *q)
{
   assert(q->D == SIMPLEX);
   SimplexParams simplexParams;
   simplexParams.deg = q->deg;
   simplexParams.dim = q->dim;
   SimplexParams *params = &simplexParams;
   BasisInterface interface = SetSimplexBasisInterface();
   Basis *simplex = BasisInit((void *)params, &interface);
   q->basis = simplex;
}


static void SetSimplexConstr(quadrature *q)
{
   assert(q->D == SIMPLEX);
   dimParamsSimplex dimSimplex = {q->dim};
   constrInterface constrinterface = SetSimplexConstrInterface();
   q->constr = constraints_init((void *)&dimSimplex, &constrinterface);
   constraints_get(q->constr);
}


static void SetCubeSimplexBasis(quadrature *q)
{
   assert(q->D == CUBESIMPLEX);
   CubeSimplexParams csParams;
   csParams.deg = q->deg;
   csParams.dims[0] = q->dims[0];
   csParams.dims[1] = q->dims[1];
   CubeSimplexParams *params = &csParams;
   BasisInterface interface = SetCubeSimplexBasisInterface();
   Basis *cubesimplex = BasisInit((void *)params, &interface);
   q->basis = cubesimplex;
}


static void SetCubeSimplexConstr(quadrature *q)
{
   assert(q->D == CUBESIMPLEX);
   dimParamsCubeSimplex dimsCubeSimplex;
   dimsCubeSimplex.dims[0] = q->dims[0];
   dimsCubeSimplex.dims[1] = q->dims[1];
   constrInterface constrinterface = SetCubeSimplexConstrInterface();
   q->constr = constraints_init((void *)&dimsCubeSimplex, &constrinterface);
   constraints_get(q->constr);
}


static void SetSimplexSimplexBasis(quadrature *q)
{
   assert(q->D == SIMPLEXSIMPLEX);
   SimplexSimplexParams ssParams;
   ssParams.deg = q->deg;
   ssParams.dims[0] = q->dims[0];
   ssParams.dims[1] = q->dims[1];
   SimplexSimplexParams *params = &ssParams;
   BasisInterface interface = SetSimplexSimplexBasisInterface();
   Basis *simplexsimplex = BasisInit((void *)params, &interface);
   q->basis = simplexsimplex;
}


static void SetSimplexSimplexConstr(quadrature *q)
{
   assert(q->D == SIMPLEXSIMPLEX);
   dimParamsSimplexSimplex dimsSimplexSimplex;
   dimsSimplexSimplex.dims[0] = q->dims[0];
   dimsSimplexSimplex.dims[1] = q->dims[1];
   constrInterface constrinterface = SetSimplexSimplexConstrInterface();
   q->constr = constraints_init((void *)&dimsSimplexSimplex, &constrinterface);
   constraints_get(q->constr);
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
      default:                 return 0;
   }
}

