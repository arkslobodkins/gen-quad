#ifndef QUADRATURE_H
#define QUADRATURE_H

#include "GENERAL_QUADRATURE.h"

#ifdef __cplusplus
extern "C" {
#endif

bool QuadInDomain(const_quadrature *q);
bool QuadPosWeights(const_quadrature *q);
bool QuadInConstraint(const_quadrature *q);
bool QuadInDomainElem(const_quadrature *q, int elem);
bool QuadInConstraintElem(const_quadrature *quad, int elem);
double QuadTestIntegral(const_quadrature *q);

void SetIntervalFuncs(quadrature *q);
void SetCubeFuncs(quadrature *q);
void SetSimplexFuncs(quadrature *q);
void SetCubeSimplexFuncs(quadrature *q);
void SetSimplexSimplexFuncs(quadrature *q);
void SetCubeSimplexSimplexFuncs(quadrature *q);

quadrature * quadrature_init(int n, int dim, int *dims, int p, DOMAIN_TYPE D);
void quad_set_funcs_and_constr(quadrature *q);
void quadrature_realloc(int n, int dim, int *dims, int p, quadrature *q);
void quadrature_reinit(int n, quadrature *q);
void quadrature_assign(const quadrature *quad1, quadrature *quad2);
void quadrature_remove_element(int index, quadrature *quad);
void quadrature_to_vector(const quadrature quad, Vector v);
void quadrature_get_elem(const_quadrature *q, int i, Vector v);
void vector_to_quadrature(const Vector v, quadrature quad);
void quadrature_free(quadrature *quad);

#ifdef __cplusplus
}
#endif

#endif
