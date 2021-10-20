#ifndef GET_JACOBIAN_H
#define GET_JACOBIAN_H

#include "GENERAL_QUADRATURE.h"
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

void GetJacobian(const int_fast8_t *basis, const quadrature quad, double *jacobian);

#ifdef __cplusplus
}
#endif

#endif
