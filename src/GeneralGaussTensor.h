#ifndef GENERALIZED_TENSOR_H
#define GENERALIZED_TENSOR_H

#include "Quadrature.h"

#ifdef __cplusplus
extern "C" {
#endif

void NodesTensor2D(const_quadrature *quad1, const_quadrature *quad2, quadrature *quad_new);
void WeightsTensor2D(const_quadrature *quad1, const_quadrature *quad2, quadrature *quad_new);
void GeneralizedNodesTensor(const_quadrature *quad1, quadrature *quad_gen);
void GeneralizedWeightsTensor(const_quadrature *quad_1, quadrature *quad_gen);

#ifdef __cplusplus
}
#endif

#endif
