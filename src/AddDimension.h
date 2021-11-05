#ifndef ADD_DIMENSION_H
#define ADD_DIMENSION_H

#include "GENERAL_QUADRATURE.h"

#ifdef __cplusplus
extern "C" {
#endif

void AddLineFirst(const quadrature *q_gauss, const quadrature *quad_prev, quadrature *quad_new);
void AddLineSimplex(const quadrature *q_gauss, const quadrature *quad_prev, quadrature *quad_new);
void GeneralDuffy(quadrature *quad);
void MixedTensor(const quadrature *q1, const quadrature *q2, quadrature *q_pr);

#ifdef __cplusplus
}
#endif

#endif
