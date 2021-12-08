#ifndef LEAST_SQUARES_NEWTONPLASMA_H
#define LEAST_SQUARES_NEWTONPLASMA_H

#include "GENERAL_QUADRATURE.h"

#ifdef __cplusplus
extern "C" {
#endif

#define SOL_FOUND true
#define SOL_NOT_FOUND false

bool LeastSquaresNewton(const bool_enum CONSTR_OPT, quadrature *q_orig, int *its);

#ifdef __cplusplus
}
#endif

#endif


