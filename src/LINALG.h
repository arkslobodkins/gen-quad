#ifndef LIN_ALG_H
#define LIN_ALG_H

#include "GENERAL_QUADRATURE.h"

#ifdef __cplusplus
extern "C" {
#endif

void dgemm_(char *TRANSA, char *TRANSB, int *M, int *N, int *K, double *ALPHA,
            double *A, int *LDA, double *B, int *LDB, double *BETA, double *C, int *LDC);


// Wrappers for LAPACK and PLASMA. All LAPACK routines are
// performed in serial, whereas PLASMA routines use OMP with
// omp_get_max_threads(), i.e. maximum number of available threads.

int DGEMM_LAPACK(CMatrix A, CMatrix B, CMatrix C);
int DGEQR2_LAPACK(CMatrix A, Vector TAU);
int DGEQRF_LAPACK(CMatrix A, Vector TAU);
int DORMQR_LAPACK(char SIDE, char TRANS, Vector TAU, CMatrix Q, CMatrix A);
int DORGQR_LAPACK(CMatrix Q, Vector TAU);
int DGESVD_LAPACK(CMatrix A, CMatrix VT);
int DGELS_LAPACK(CMatrix A, Vector RHS_TO_X);

#ifdef _OPENMP
int DGELS_PLASMA(CMatrix A, Vector RHS_TO_X);
int DGEMM_PLASMA(CMatrix A, CMatrix B, CMatrix C);
#endif

void Transpose(int M, int N, const double *A, double *B);

#ifdef __cplusplus
}
#endif

#endif
