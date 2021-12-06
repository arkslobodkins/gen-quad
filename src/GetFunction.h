#ifndef GET_FUNCTION_H
#define GET_FUNCTION_H

#include "GENERAL_QUADRATURE.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _OPENMP
void GetFunctionOmp(quadrature *q, Vector f);
#endif
void GetFunction(quadrature *q, Vector f);

#ifdef __cplusplus
}
#endif

#endif
