#ifndef SET_PARAMS_H
#define SET_PARAMS_H

#include "GENERAL_QUADRATURE.h"

#ifdef __cplusplus
extern "C" {
#endif

void SetParams(int dim, int num_dims, int *dims, int deg, quadParams *params);

#ifdef __cplusplus
}
#endif

#endif
