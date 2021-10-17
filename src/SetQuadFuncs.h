#ifndef SET_QUAD_FUNCS_H
#define SET_QUAD_FUNCS_H

#include "GENERAL_QUADRATURE.h"

#ifdef __cplusplus
extern "C" {
#endif

void SetIntervalFuncs(quadrature *q);
void SetCubeFuncs(quadrature *q);
void SetSimplexFuncs(quadrature *q);
void SetCubeSimplexFuncs(quadrature *q);
void SetSimplexSimplexFuncs(quadrature *q);
void SetCubeSimplexSimplexFuncs(quadrature *q);

#ifdef __cplusplus
}
#endif

#endif
