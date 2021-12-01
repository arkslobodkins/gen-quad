#ifndef BASIS_FUNCTIONS_H
#define BASIS_FUNCTIONS_H

#include "GENERAL_QUADRATURE.h"

#ifdef __cplusplus
extern "C" {
#endif

void LegendrePoly(int order, double x, double *p, double *dp);
void JacobiPoly(int order, double x, double alpha, double beta, double *p);
void JacobiPolyPrime(int order, double x, double alpha, double beta, double *dp);

void CubeBasisFuncs(int *dims, int deg, const INT_8 *basis_id, const double *x, double *phi);
void CubeBasisFuncsDer(int *dims, int deg, const INT_8 *basis_id, const double *x, double *phiPrime);

void BasisSimplex(int *dims, int deg, const INT_8 *basis_id, const double *x, double *phi);
void BasisPrimeSimplex(int *dims, int deg, const INT_8 *basis_id, const double *x, double *phiPrime);

void BasisCubeSimplex(int *dims, int deg, const INT_8 *basis_id, const double *x, double *phi);
void BasisPrimeCubeSimplex(int *dims, int deg, const INT_8 *basis_id, const double *x, double *phiPrime);

void BasisSimplexSimplex(int *dims, int deg, const INT_8 *basis_id, const double *x, double *phi);
void BasisPrimeSimplexSimplex(int *dims, int deg, const INT_8 *basis_id, const double *x, double *phiPrime);

void BasisMonomial(int dim, int deg, const INT_8 *basis_id, const double *x, double *phi);

void orthogonal_simplex_basis_test(int deg, int dim);

#ifdef __cplusplus
}
#endif

#endif

