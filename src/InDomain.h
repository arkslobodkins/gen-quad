#ifndef IN_DOMAIN_TYPE_H
#define IN_DOMAIN_TYPE_H

#include "GENERAL_QUADRATURE.h"
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

bool InCubeElem(const_quadrature *q, int elem);
bool InSimplexElem(const_quadrature *q, int elem);
bool InCubeSimplexElem(const_quadrature *q, int elem);
bool InSimplexSimplexElem(const_quadrature *q, int elem);
bool InCubeSimplexSimplexElem(const_quadrature *q, int elem);

bool InCube(const_quadrature *q);
bool InSimplex(const_quadrature *q);
bool InCubeSimplex(const_quadrature *q);
bool InSimplexSimplex(const_quadrature *q);
bool InCubeSimplexSimplex(const_quadrature *q);

#ifdef __cplusplus
}
#endif

#endif
