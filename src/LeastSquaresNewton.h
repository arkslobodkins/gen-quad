#ifndef LEAST_SQUARES_NEWTON_H
#define LEAST_SQUARES_NEWTON_H

#include "GENERAL_QUADRATURE.h"
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

#define SOL_FOUND true
#define SOL_NOT_FOUND false

bool LeastSquaresNewton(const bool_enum FLAG_CONSTR, const INT_8 *basis, quadrature *q_orig, int *its);

#ifdef __cplusplus
}
#endif

#endif


