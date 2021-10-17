#ifndef JACOBI_POLY_H
#define JACOBI_POLY_H

#ifdef __cplusplus
extern "C" {
#endif

void JacobiPoly(int order, double x, double alpha, double beta, double *p);
void JacobiPolyPrime(int order, double x, double alpha, double beta, double *dp);

#ifdef __cplusplus
}
#endif

#endif
