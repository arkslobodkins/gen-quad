#ifndef BASIS_INDICES_H
#define BASIS_INDICES_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

void BasisIndices(int deg, int dim, int_fast8_t *f);
int BasisSize(int deg, int dim);

#ifdef __cplusplus
}
#endif

#endif
