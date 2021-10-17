#ifndef GAUSS_TENSOR_H
#define GAUSS_TENSOR_H

#include "Quadrature.h"

#ifdef __cplusplus
extern "C" {
#endif

void NodesTensor(const_quadrature *quad1, const_quadrature *quad2, quadrature *quad_new);
void WeightsTensor(const_quadrature *quad1, const_quadrature *quad2, quadrature *quad_new);

#ifdef __cplusplus
}
#endif

#endif
