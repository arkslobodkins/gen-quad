/* Arkadijs Slobodkins
 * SMU Mathematics
 * Copyright 2022
 */


#ifndef LIN_ALG_H
#define LIN_ALG_H

#include "GENERAL_QUADRATURE.h"

#ifdef __cplusplus
extern "C" {
#endif

// Wrappers for MKL LAPACK and PLASMA. All LAPACK routines are
// performed in serial, whereas PLASMA routines use OpenMP with
// omp_get_max_threads(), i.e. maximum number of available threads.

void DGEMM_LAPACK(CMatrix A, CMatrix B, CMatrix C);
int DGEQR2_LAPACK(CMatrix A, Vector TAU);
int DGEQRF_LAPACK(CMatrix A, Vector TAU);
int DORMQR_LAPACK(char SIDE, char TRANS, Vector TAU, CMatrix Q, CMatrix A);
int DORGQR_LAPACK(CMatrix Q, Vector TAU);
int DGESVD_LAPACK(CMatrix A, CMatrix VT);
int DGELS_LAPACK(CMatrix A, Vector b, Vector x);

#ifdef _OPENMP
int DGEMM_PLASMA(CMatrix A, CMatrix B, CMatrix C);
int DGELS_PLASMA(CMatrix A, Vector b, Vector x);
#endif

#ifdef __cplusplus
}
#endif

#endif
