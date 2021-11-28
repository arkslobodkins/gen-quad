#ifndef LIN_ALG_H
#define LIN_ALG_H

#include "GENERAL_QUADRATURE.h"

#ifdef __cplusplus
extern "C" {
#endif

void dgemm_(char *TRANSA, char *TRANSB, int *M, int *N, int *K, double *ALPHA,
            double *A, int *LDA, double *B, int *LDB, double *BETA, double *C, int *LDC);

void dgels_(char *TRANS, int *M, int *N, int *NRHS,
            double *A, int *LDA, double *B, int *LDB,
            double *WORK, int *LWORK, int *INFO);

void dgeqr2_(int * M, int *N, double *A, int *LDA,
             double *TAU, double *WORK, int *INFO);

void dgeqrf_(int * M, int *N, double *A, int *LDA,
             double *TAU, double *WORK, int *LWORK, int *INFO);

void dormqr_(char *SIDE, char *TRANS, int *M, int *N,
             int *size_TAU, double *Q, int *LDQ, double *TAU,
             double *C, int *LDC,
             double *WORK, int *LWORK, int *INFO);

// Wrappers for LAPACK and PLASMA. All LAPACK routines are
// performed in serial, whereas PLASMA routines use OMP with
// omp_get_max_threads(), i.e. maximum number of available threads.

void MATMUL_LAPACK(int M, int K, int N, double *MAT1, double *MAT2, double *MAT3);

int DGEQR2_LAPACK(CMatrix A, Vector TAU);
int DGEQRF_LAPACK(CMatrix A, Vector TAU);
int DORMQR_LAPACK(char SIDE, char TRANS, Vector TAU, CMatrix Q, CMatrix A);

int DGELS_LAPACK(CMatrix A, Vector RHS_TO_X);
#ifdef _OPENMP
int DGELS_PLASMA(CMatrix A, Vector RHS_TO_X);
#endif

void Transpose(int M, int N, const double *A, double *B); // hand-written

#ifdef __cplusplus
}
#endif

#endif
