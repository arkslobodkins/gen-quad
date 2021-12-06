#ifndef QUADRATURE_H
#define QUADRATURE_H

#include "GENERAL_QUADRATURE.h"
#include "Basis.h"

#ifdef __cplusplus
extern "C" {
#endif


quadrature * quadrature_init_basic(int n, int dim, int *dims, int deg, DOMAIN_TYPE D);
quadrature * quadrature_init_full(int n, int dim, int *dims, int deg, DOMAIN_TYPE D);
void quadrature_realloc_array(int n, quadrature *q);
void quadrature_reinit_simple(int n, int dim, int *dims, int deg, quadrature *q);
void quadrature_shrink_array(int n, quadrature *q);
quadrature *quadrature_make_full_copy(const quadrature *q);
quadrature *quadrature_without_element(quadrature *q, int i);
void quadrature_assign(const quadrature *quad1, quadrature *quad2);
void quadrature_remove_element(int index, quadrature *quad);
void quadrature_to_vector(const quadrature quad, Vector v);
void quadrature_get_elem(const quadrature *q, int i, Vector v);
void vector_to_quadrature(const Vector v, quadrature quad);
void quadrature_free(quadrature *quad);

#ifdef _OPENMP
void QuadAllocBasisOmp(quadrature *q, int num_threads);
void QuadFreeBasisOmp(quadrature *q, int num_threads);
#endif
bool QuadInDomain(const quadrature *q);
bool QuadInDomainElem(const quadrature *q, int elem);
bool QuadPosWeights(const quadrature *q);
bool QuadInConstraint(const quadrature *q);
bool QuadInConstraintElem(const quadrature *quad, int elem);

bool QuadOnTheBoundary(const quadrature *q, int elem);
bool QuadEqnOnTheBoundary(const quadrature *q, int elem, int eqn);

double QuadTestIntegral(const quadrature *q, BASIS_TYPE btype);
double QuadTestIntegralExp(const quadrature *q);

#ifdef __cplusplus
}
#endif

#endif
