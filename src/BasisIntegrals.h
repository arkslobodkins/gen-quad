#ifndef INTEGRALS_OF_BASIS_FUNCTIONS_H
#define INTEGRALS_OF_BASIS_FUNCTIONS_H

#include "GENERAL_QUADRATURE.h"

#ifdef __cplusplus
extern "C" {
#endif

void IntegralsCube(const quadParams *params, double *integrals);
void IntegralsSimplex(const quadParams *params, double *integrals);
void IntegralsCubeSimplex(const quadParams *params, double *integrals);
void IntegralsSimplexSimplex(const quadParams *params, double *integrals);
void IntegralsCubeSimplexSimplex(const quadParams *params, double *integrals);

#ifdef __cplusplus
}
#endif

#endif

