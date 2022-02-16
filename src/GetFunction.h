#ifndef GET_FUNCTION_H
#define GET_FUNCTION_H

#include "GENERAL_QUADRATURE.h"
#include "Quadrature.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _OPENMP
void GetFunctionOmp(quadrature *q, Vector f);
#endif
void GetFunction(quadrature *q, Vector f);
void TestResidual(quadrature *q, const char* str);

#ifdef __cplusplus
}
#endif

#endif
