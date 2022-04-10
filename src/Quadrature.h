/* Arkadijs Slobodkins
 * SMU Mathematics
 * Copyright 2022
 */


#ifndef QUADRATURE_H
#define QUADRATURE_H

#include "Basis.h"

#ifdef __cplusplus
extern "C" {
#endif


typedef struct quadrature quadrature;
typedef void(*SetBasisAndConstr)(quadrature *q);
typedef void(*SetBasis)(quadrature *q);
typedef void(*SetConstr)(quadrature *q);
typedef double(*ExpIntegralExact)(const quadrature *q);
typedef void(*FreePtr)(quadrature *quad);


struct quadrature
{
   DOMAIN_TYPE D;
   int dim;
   int num_dims;
   int *dims;
   int deg;

   int num_nodes;
   double *w;
   double *x;
   Vector z;

   int isFullyInitialized;
   SetBasisAndConstr setBasisAndConstr;
   SetBasis setBasis;
   SetConstr setConstr;

   struct constraints *constr;
   struct Basis *basis;
   struct Basis **basisOmp;
   ExpIntegralExact expIntegralExact;

   FreePtr free_ptr;
};

static inline double qweight(const quadrature *q, int i)
{
   assert(i >= 0 && i <= q->num_nodes-1);
   return q->w[i];
}

static inline double* qnode(const quadrature *q, int i)
{
   assert(i >= 0 && i <= q->num_nodes-1);
   return &q->x[q->dim*i];
}


quadrature* quadrature_init_basic(int n, int dim, int *dims, int deg, DOMAIN_TYPE D);
quadrature* quadrature_init_full(int n, int dim, int *dims, int deg, DOMAIN_TYPE D);
quadrature* quadrature_gauss_legendre(int deg);
quadrature* quadrature_gauss_jacobi(int deg, double alpha, double beta);
quadrature* quadrature_full_cube_tensor(int deg, int dim);
quadrature* quadrature_full_simplex_tensor(int deg, int dim);
void quadrature_reinit_basic(int n, int dim, int *dims, int deg, quadrature *q);
void quadrature_realloc_array(int n, quadrature *q);

// Routine that reallocates q to a quadrature with n nodes.
// First n nodes are copied to the new quadrature.
// Before the routine is called, q should contain at least n nodes or more.
void quadrature_shrink_array(int n, quadrature *q);

quadrature* quadrature_make_full_copy(const quadrature *q);
quadrature* quadrature_without_element(const quadrature *q, int i);
void quadrature_assign(const quadrature *quad1, quadrature *quad2);
void quadrature_assign_resize(const quadrature *q1, quadrature *q2);
void quadrature_remove_element(int index, quadrature *quad);
void quadrature_to_vector(const quadrature quad, Vector v);
void quadrature_get_elem(const quadrature *q, int i, Vector v);
void vector_to_quadrature(const Vector v, quadrature quad);
void quadrature_fill_random(quadrature *q);
void quadrature_free(quadrature *quad);

#ifdef _OPENMP
void QuadAllocBasisOmp(quadrature *q, int num_threads);
void QuadFreeBasisOmp(quadrature *q, int num_threads);
#endif
bool QuadInDomainElem(const quadrature *q, int elem);
bool QuadInDomainElemEps(const quadrature *q, int elem);
bool QuadInDomain(const quadrature *q);
bool QuadInDomainEps(const quadrature *q);
bool QuadInDomainEqnElemEps(const quadrature *q, int elem, int eqn);

bool QuadPosWeights(const quadrature *q);
bool QuadPosWeightsEps(const quadrature *q);

bool QuadInConstraintElem(const quadrature *quad, int elem);
bool QuadInConstraint(const quadrature *q);
bool QuadInConstraintEps(const quadrature *q);

bool QuadOnTheBoundary(const quadrature *q, int elem);
bool QuadEqnOnTheBoundary(const quadrature *q, int elem, int eqn);
double QuadDistFromTheBoundaryElem(const quadrature *q, int elem);
double QuadMinDistFromTheBoundary(const quadrature *q);
double QuadDistFromTheBoundaryTwoNorm(const quadrature *q);

double QuadDistAlphaElem(const quadrature *q, int elem);
VMin QuadMinAlpha(const quadrature *q);

bool QuadIsEqual(const quadrature *a, const quadrature *b);
double QuadTestIntegral(const quadrature *q, BASIS_TYPE btype);
double QuadTestIntegralExp(const quadrature *q);

#ifdef __cplusplus
}
#endif

#endif
