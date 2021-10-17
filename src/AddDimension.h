#ifndef ADD_DIMENSION_H
#define ADD_DIMENSION_H

#include "GENERAL_QUADRATURE.h"

#ifdef __cplusplus
extern "C" {
#endif

void AddLineFirst(const_quadrature *q_gauss, const_quadrature *quad_prev, quadrature *quad_new);
void AddLineSimplex(const_quadrature *q_gauss, const_quadrature *quad_prev, quadrature *quad_new);
void GeneralDuffy(quadrature *quad);

#ifdef __cplusplus
}
#endif

#endif
