#ifndef GET_JACOBIAN_H
#define GET_JACOBIAN_H

#include "GENERAL_QUADRATURE.h"

#ifdef __cplusplus
extern "C" {
#endif

void GetJacobian(const INT_8 *basisIndices, quadrature *q, CMatrix);

#ifdef __cplusplus
}
#endif

#endif
