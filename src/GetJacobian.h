#ifndef GET_JACOBIAN_H
#define GET_JACOBIAN_H

#include "GENERAL_QUADRATURE.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _OPENMP
void GetJacobianOmp(quadrature *q, CMatrix);
#endif
void GetJacobian(quadrature *q, CMatrix);

#ifdef __cplusplus
}
#endif

#endif
