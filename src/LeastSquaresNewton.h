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

#define CONSTRAINT_SUCCESS true
#define CONSTRAINT_FAILURE false

bool LeastSquaresNewton(const bool_enum FLAG_CONSTR, const int_fast8_t *basis, quadrature *q_orig, int *its);

ConstrNodeData constrain_vector(const Matrix A, const Vector b, const_quadrature *q_prev, const_quadrature *q_next);
#ifdef __cplusplus
}
#endif

#endif


