#ifndef COMPUTE_DOMAIN_TYPE_H
#define COMPUTE_DOMAIN_TYPE_H

#include "GENERAL_QUADRATURE.h"

#ifdef __cplusplus
extern "C" {
#endif

void ComputeCubeFull(int deg, int dim);

void ComputeInterval(int degree);
void ComputeCube(int degree, int dim);
void ComputeSimplex(int degree, int dim);
void ComputeCubeSimplex(int degree, int dim1, int dim2);
void ComputeSimplexSimplex(int degree, int dim1, int dim2);

#ifdef __cplusplus
}
#endif

#endif
