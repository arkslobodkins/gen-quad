/* Arkadijs Slobodkins
 * SMU Mathematics
 * Copyright 2022
 */


#ifndef GET_JACOBIAN_H
#define GET_JACOBIAN_H

#include "GENERAL_QUADRATURE.h"
#include "Quadrature.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _OPENMP
void GetJacobianOmp(quadrature *q, CMatrix);
void GetFunctionAndJacobianOmp(quadrature *q, Vector f, CMatrix JACOBIAN);
#endif

// Computes an returns Jacobian ∈  R[num_funcs, (dim+1)*num_nodes]
// of a vector function F = P(X) * W - In,
// P(X) ∈  R[numFuncs, num_nodes],  W ∈  R[num_nodes], In ∈ R[numFuncs],
// where each row of P(X) corresponds to a distinct polynomial of orthogonal basis,
// each column of P(X) and entries in W correspond to ith node and weight of the quadrature,
// and In are integrals of the polynomial basis over the domain specified by quadrature.
void GetJacobian(quadrature *q, CMatrix);

void GetFunctionAndJacobian(quadrature *q, Vector f, CMatrix JACOBIAN);
void GetBasis(quadrature *q, CMatrix BasisMatrix);

#ifdef __cplusplus
}
#endif

#endif
