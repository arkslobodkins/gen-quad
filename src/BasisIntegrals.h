#ifndef INTEGRALS_OF_BASIS_FUNCTIONS_H
#define INTEGRALS_OF_BASIS_FUNCTIONS_H

#include "GENERAL_QUADRATURE.h"

#ifdef __cplusplus
extern "C" {
#endif

void BasisIntegralsCube(int *dims, int deg, double *integrals);
void BasisIntegralsSimplex(int *dims, int deg, double *integrals);
void BasisIntegralsCubeSimplex(int *dims, int deg, double *integrals);
void BasisIntegralsSimplexSimplex(int *dims, int deg, double *integrals);
void BasisIntegralsCubeSimplexSimplex(int *dims, int deg, double *integrals);

#ifdef __cplusplus
}
#endif

#endif

