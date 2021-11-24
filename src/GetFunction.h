#ifndef GET_FUNCTION_H
#define GET_FUNCTION_H

#include "GENERAL_QUADRATURE.h"

#ifdef __cplusplus
extern "C" {
#endif

void GetFunction(const INT_8 *basisIndices, quadrature *q, Vector f);

#ifdef __cplusplus
}
#endif

#endif
