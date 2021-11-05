#ifndef LEAST_SQUARES_NEWTON_H
#define LEAST_SQUARES_NEWTON_H

#include "GENERAL_QUADRATURE.h"
#include <stdbool.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define SOL_FOUND true
#define SOL_NOT_FOUND false

#define CONSTR_SUCCESS  0
#define CONSTR_FAILURE -1

bool LeastSquaresNewton(const bool_enum FLAG_CONSTR, const int_fast8_t *basis, quadrature *q_orig, int *its);

#ifdef CONSTR_OPT
ConstrVectData ConstrVectDataInit();
int ConstrainedOptimization(const quadrature *q_prev, quadrature *q_next, ConstrVectData *cVecData);
#endif

#ifdef __cplusplus
}
#endif

#endif


