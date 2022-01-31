#ifndef GET_JACOBIAN_H
#define GET_JACOBIAN_H

#include "GENERAL_QUADRATURE.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _OPENMP
void GetJacobianOmp(quadrature *q, CMatrix);
void GetFunctionAndJacobianOmp(quadrature *q, Vector f, CMatrix JACOBIAN);
#endif
void GetJacobian(quadrature *q, CMatrix);
void GetFunctionAndJacobian(quadrature *q, Vector f, CMatrix JACOBIAN);
void GetBasis(quadrature *q, CMatrix BasisMatrix);

#ifdef __cplusplus
}
#endif

#endif
