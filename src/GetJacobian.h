#ifndef GET_JACOBIAN_H
#define GET_JACOBIAN_H

#include "GENERAL_QUADRATURE.h"
#include "Quadrature.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _OPENMP
// Computes an returns Jacobian of size num_funcs x (dim+1)*num_nodes
// of a vector X * w, where each row of X corresponds to
// a distinct polynomial of orthogonal basis, and each
// column of X and entries in w correspond to ith node and weight
// of the quadrature.
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
