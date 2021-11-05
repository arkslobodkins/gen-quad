#ifndef LIN_ALG_H
#define LIN_ALG_H

#include "Matrix.h"

#ifdef __cplusplus
extern "C" {
#endif

void dgemm_(char *, char *, int *, int *, int *, double *, double *, int *, double *, int *, double *, double *, int *);
void dgeqr2_(int *, int *, double *, int *, double *, double *, int *);
void dormqr_(char *, char *, int *, int *, int *, double *, int *, double *, double *, int *, double *, int *, int *);
void dgels_(char *, int *, int *, int *, double *, int *, double*, int *, double *, int *, int *);

void MAT_MUL(int M, int K, int N, double *MAT1, double *MAT2, double *MAT3);
int DORMQR_M(char SIDE, char TRANS, double *REFL_Q, CMatrix Q, CMatrix A);
int DORMQR_V(char SIDE, char TRANS, double *REFL_Q, CMatrix Q, Vector A);
void Transpose(int M, int N, const double *A, double *B);

#ifdef __cplusplus
}
#endif

#endif
