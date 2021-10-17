#ifndef PHI_H
#define PHI_H

#include "GENERAL_QUADRATURE.h"
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

void PhiCube(const int_fast8_t *basis_id, const double *x,
             const quadParams *params, double *phi);

void PhiPrimeCube(const int_fast8_t *basis_id, const double *x,
                  const quadParams *params, double *phiPrime);

void PhiSimplex(const int_fast8_t *basis_id, const double *x,
                const quadParams *params, double *phi);

void PhiPrimeSimplex(const int_fast8_t *basis_id, const double *x,
                     const quadParams *params, double *phiPrime);

void PhiCubeSimplex(const int_fast8_t *basis_id, const double *x,
                    const quadParams *params, double *phi);

void PhiPrimeCubeSimplex(const int_fast8_t *basis_id, const double *x,
                         const quadParams *params, double *phiPrime);

void PhiSimplexSimplex(const int_fast8_t *basis_id, const double *x,
                       const quadParams *params, double *phi);

void PhiPrimeSimplexSimplex(const int_fast8_t *basis_id, const double *x,
                            const quadParams *params, double *phiPrime);

void PhiCubeSimplexSimplex(const int_fast8_t *basis_id, const double *x,
                           const quadParams *params, double *phi);

void PhiPrimeCubeSimplexSimplex(const int_fast8_t *basis_id, const double *x,
                                const quadParams *params, double *phiPrime);

void PhiMonomial(const int_fast8_t *basis_id, const double *x,
                 const quadParams *params, double *phi);

void PhiPrimeMonomial(const int_fast8_t *basis_id, const double *x,
                      const quadParams *params, double *phiPrime);

double  orthogonal_simplex_basis_test(const int_fast8_t *basis_id, const quadParams *params);

#ifdef __cplusplus
}
#endif

#endif

