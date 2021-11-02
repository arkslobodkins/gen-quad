#ifndef GET_FUNCTION_H
#define GET_FUNCTION_H

#include "GENERAL_QUADRATURE.h"
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

void GetFunction(const int_fast8_t *basis, quadrature *q, Vector f);

#ifdef __cplusplus
}
#endif

#endif
