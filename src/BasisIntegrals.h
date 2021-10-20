#ifndef INTEGRALS_OF_BASIS_FUNCTIONS_H
#define INTEGRALS_OF_BASIS_FUNCTIONS_H

#include "GENERAL_QUADRATURE.h"

#ifdef __cplusplus
extern "C" {
#endif

void IntegralsCube(int *dims, int deg, double *integrals);
void IntegralsSimplex(int *dims, int deg, double *integrals);
void IntegralsCubeSimplex(int *dims, int deg, double *integrals);
void IntegralsSimplexSimplex(int *dims, int deg, double *integrals);
void IntegralsCubeSimplexSimplex(int *dims, int deg, double *integrals);

#ifdef __cplusplus
}
#endif

#endif

