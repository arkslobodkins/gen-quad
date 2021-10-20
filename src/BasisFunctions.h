#ifndef BASIS_FUNCTIONS_H
#define BASIS_FUNCTIONS_H

#include "GENERAL_QUADRATURE.h"
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

void PhiCube(int *dims, int deg, const int_fast8_t *basis_id, const double *x, double *phi);
void PhiPrimeCube(int *dims, int deg, const int_fast8_t *basis_id, const double *x, double *phiPrime);

void PhiSimplex(int *dims, int deg, const int_fast8_t *basis_id, const double *x, double *phi);
void PhiPrimeSimplex(int *dims, int deg, const int_fast8_t *basis_id, const double *x, double *phiPrime);

void PhiCubeSimplex(int *dims, int deg, const int_fast8_t *basis_id, const double *x, double *phi);
void PhiPrimeCubeSimplex(int *dims, int deg, const int_fast8_t *basis_id, const double *x, double *phiPrime);

void PhiSimplexSimplex(int *dims, int deg, const int_fast8_t *basis_id, const double *x, double *phi);
void PhiPrimeSimplexSimplex(int *dims, int deg, const int_fast8_t *basis_id, const double *x, double *phiPrime);

void PhiCubeSimplexSimplex(int *dims, int deg, const int_fast8_t *basis_id, const double *x, double *phi);
void PhiPrimeCubeSimplexSimplex(int *dims, int deg, const int_fast8_t *basis_id, const double *x, double *phiPrime);

double  orthogonal_simplex_basis_test(int *dims, int deg, const int_fast8_t *basis_id);

#ifdef __cplusplus
}
#endif

#endif

