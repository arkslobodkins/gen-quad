#ifndef ADD_DIMENSION_H
#define ADD_DIMENSION_H

#include "GENERAL_QUADRATURE.h"

#ifdef __cplusplus
extern "C" {
#endif

void AddLineFirst(const quadrature *q1D, const quadrature *quad_prev, quadrature *quad_new);
void AddLineSimplex(const quadrature *q1D, const quadrature *quad_prev, quadrature *quad_new);
void GeneralDuffy(quadrature *q);
void MixedTensor(const quadrature *q1, const quadrature *q2, quadrature *q_tp);

#ifdef __cplusplus
}
#endif

#endif
