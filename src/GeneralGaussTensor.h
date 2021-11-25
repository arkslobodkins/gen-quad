#ifndef GENERALIZED_TENSOR_H
#define GENERALIZED_TENSOR_H

#include "Quadrature.h"

#ifdef __cplusplus
extern "C" {
#endif

void NodesTensor2D(const quadrature *quad1, const quadrature *quad2, quadrature *quad_new);
void WeightsTensor2D(const quadrature *quad1, const quadrature *quad2, quadrature *quad_new);
void GeneralizedNodesTensor(const quadrature *quad1, quadrature *quad_gen);
void GeneralizedWeightsTensor(const quadrature *quad_1, quadrature *quad_gen);
void TestGeneralizedTensor();

#ifdef __cplusplus
}
#endif

#endif
