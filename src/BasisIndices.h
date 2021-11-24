#ifndef BASIS_INDICES_H
#define BASIS_INDICES_H

#include "GENERAL_QUADRATURE.h"

#ifdef __cplusplus
extern "C" {
#endif

void BasisIndices(int deg, int dim, INT_8 *f);
int BasisSize(int deg, int dim);

#ifdef __cplusplus
}
#endif

#endif
