#ifndef LIN_ALG_H
#define LIN_ALG_H

#include "GENERAL_QUADRATURE.h"

#ifdef __cplusplus
extern "C" {
#endif

// LAPACK prototypes

void dorgqr_(int *M, int *N, int *K, double *Q, int *LDA,
             double *TAU, double *WORK, int *LWORK, int *INFO);
void dgemm_(char *TRANSA, char *TRANSB, int *M, int *N, int *K, double *ALPHA,
            double *A, int *LDA, double *B, int *LDB, double *BETA, double *C, int *LDC);

void dgels_(char *TRANS, int *M, int *N, int *NRHS,
            double *A, int *LDA, double *B, int *LDB,
            double *WORK, int *LWORK, int *INFO);

void dgeqr2_(int *M, int *N, double *A, int *LDA,
             double *TAU, double *WORK, int *INFO);

void dgeqrf_(int *M, int *N, double *A, int *LDA,
             double *TAU, double *WORK, int *LWORK, int *INFO);

void dormqr_(char *SIDE, char *TRANS, int *M, int *N,
             int *SIZE_TAU, double *Q, int *LDQ, double *TAU,
             double *C, int *LDC,
             double *WORK, int *LWORK, int *INFO);

void dgesvd_(char *JOBU, char *JOBVT,
		       int *M, int *N, double *A, int *LDA,
             double *S, double *U, int *LDU, double *VT, int *LDVT,
             double *WORK, int *LWORK, int *INFO);


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
void Transpose(int M, int N, const double *A, double *B);
#endif

#ifdef __cplusplus
}
#endif

#endif
