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


bool LeastSquaresNewton(const BOOLEAN FLAG_CONSTR, const int_fast8_t *basis,
                        const _DomainFuncs funcs, const constraints *cons,
                        int *its, quadrature *q_orig);

#ifdef __cplusplus
}
#endif

#endif


