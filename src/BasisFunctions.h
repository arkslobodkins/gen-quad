#ifndef BASIS_FUNCTIONS_H
#define BASIS_FUNCTIONS_H

#include "GENERAL_QUADRATURE.h"
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

void BasisCube(int *dims, int deg, const INT_8 *basis_id, const double *x, double *phi);
void BasisPrimeCube(int *dims, int deg, const INT_8 *basis_id, const double *x, double *phiPrime);

void BasisSimplex(int *dims, int deg, const INT_8 *basis_id, const double *x, double *phi);
void BasisPrimeSimplex(int *dims, int deg, const INT_8 *basis_id, const double *x, double *phiPrime);

void BasisCubeSimplex(int *dims, int deg, const INT_8 *basis_id, const double *x, double *phi);
void BasisPrimeCubeSimplex(int *dims, int deg, const INT_8 *basis_id, const double *x, double *phiPrime);

void BasisSimplexSimplex(int *dims, int deg, const INT_8 *basis_id, const double *x, double *phi);
void BasisPrimeSimplexSimplex(int *dims, int deg, const INT_8 *basis_id, const double *x, double *phiPrime);

void BasisCubeSimplexSimplex(int *dims, int deg, const INT_8 *basis_id, const double *x, double *phi);
void BasisPrimeCubeSimplexSimplex(int *dims, int deg, const INT_8 *basis_id, const double *x, double *phiPrime);

void BasisMonomial(int dim, int deg, const INT_8 *basis_id, const double *x, double *phi);

double  orthogonal_simplex_basis_test(int *dims, int deg, const INT_8 *basis_id);

#ifdef __cplusplus
}
#endif

#endif

