/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from core_blas/core_zlantr.c, normal z -> s, Thu Mar 10 18:58:42 2022
 *
 **/

#include "core_blas.h"
#include "plasma_types.h"
#include "plasma_internal.h"
#include "core_lapack.h"

#include <math.h>

/******************************************************************************/
void core_slantr(plasma_enum_t norm, plasma_enum_t uplo, plasma_enum_t diag,
                 int m, int n,
                 const float *A, int lda,
                 float *work, float *value)
{
    // Due to a bug in LAPACKE < 3.6.1, this function always returns zero.
    // *value = LAPACKE_slantr_work(LAPACK_COL_MAJOR,
    //                              lapack_const(norm), lapack_const(uplo),
    //                              lapack_const(diag),
    //                              m, n, A, lda, work);

    // Calling LAPACK directly instead.
    char nrm = lapack_const(norm);
    char upl = lapack_const(uplo);
    char dia = lapack_const(diag);
    *value = LAPACK_slantr(&nrm, &upl, &dia, &m, &n, A, &lda, work);
}

/******************************************************************************/
void core_omp_slantr(plasma_enum_t norm, plasma_enum_t uplo, plasma_enum_t diag,
                     int m, int n,
                     const float *A, int lda,
                     float *work, float *value,
                     plasma_sequence_t *sequence, plasma_request_t *request)
{
    #pragma omp task depend(in:A[0:lda*n]) \
                     depend(out:value[0:1])
    {
        if (sequence->status == PlasmaSuccess)
            core_slantr(norm, uplo, diag, m, n, A, lda, work, value);
    }
}

/******************************************************************************/
void core_omp_slantr_aux(plasma_enum_t norm, plasma_enum_t uplo,
                         plasma_enum_t diag,
                         int m, int n,
                         const float *A, int lda,
                         float *value,
                         plasma_sequence_t *sequence, plasma_request_t *request)
{
    switch (norm) {
    case PlasmaOneNorm:
        #pragma omp task depend(in:A[0:lda*n]) \
                         depend(out:value[0:n])
        {
            if (sequence->status == PlasmaSuccess) {
                if (uplo == PlasmaUpper) {
                    if (diag == PlasmaNonUnit) {
                        for (int j = 0; j < n; j++) {
                            value[j] = fabsf(A[lda*j]);
                            for (int i = 1; i < imin(j+1, m); i++) {
                                value[j] += fabsf(A[lda*j+i]);
                            }
                        }
                    }
                    else { // PlasmaUnit
                        int j;
                        for (j = 0; j < imin(n, m); j++) {
                            value[j] = 1.0;
                            for (int i = 0; i < j; i++) {
                                value[j] += fabsf(A[lda*j+i]);
                            }
                        }
                        for (; j < n; j++) {
                            value[j] = fabsf(A[lda*j]);
                            for (int i = 1; i < m; i++) {
                                value[j] += fabsf(A[lda*j+i]);
                            }
                        }
                    }
                }
                else { // PlasmaLower
                    if (diag == PlasmaNonUnit) {
                        int j;
                        for (j = 0; j < imin(n, m); j++) {
                            value[j] = fabsf(A[lda*j+j]);
                            for (int i = j+1; i < m; i++) {
                                value[j] += fabsf(A[lda*j+i]);
                            }
                        }
                        for (; j < n; j++)
                            value[j] = 0.0;
                    }
                    else { // PlasmaUnit
                        int j;
                        for (j = 0; j < imin(n, m); j++) {
                            value[j] = 1.0;
                            for (int i = j+1; i < m; i++) {
                                value[j] += fabsf(A[lda*j+i]);
                            }
                        }
                        for (; j < n; j++)
                            value[j] = 0.0;
                    }
                }
            }
        }
        break;
    case PlasmaInfNorm:
        #pragma omp task depend(in:A[0:lda*n]) \
                         depend(out:value[0:m])
        {
            if (sequence->status == PlasmaSuccess) {
                if (uplo == PlasmaUpper) {
                    if (diag == PlasmaNonUnit) {
                        for (int i = 0; i < m; i++)
                            value[i] = 0.0;

                        for (int j = 0; j < n; j++) {
                            for (int i = 0; i < imin(j+1, m); i++) {
                                value[i] += fabsf(A[lda*j+i]);
                            }
                        }
                    }
                    else { // PlasmaUnit
                        int i;
                        for (i = 0; i < imin(m, n); i++)
                            value[i] = 1.0;

                        for (; i < m; i++)
                            value[i] = 0.0;

                        int j;
                        for (j = 0; j < imin(n, m); j++) {
                            for (i = 0; i < j; i++) {
                                value[i] += fabsf(A[lda*j+i]);
                            }
                        }
                        for (; j < n; j++) {
                            for (i = 0; i < m; i++) {
                                value[i] += fabsf(A[lda*j+i]);
                            }
                        }
                    }
                }
                else { // PlasmaLower
                    if (diag == PlasmaNonUnit) {
                        for (int i = 0; i < m; i++)
                            value[i] = 0.0;

                        for (int j = 0; j < imin(n, m); j++) {
                            for (int i = j; i < m; i++) {
                                value[i] += fabsf(A[lda*j+i]);
                            }
                        }
                    }
                    else { // PlasmaUnit
                        int i;
                        for (i = 0; i < imin(m, n); i++)
                            value[i] = 1.0;

                        for (; i < m; i++)
                            value[i] = 0.0;

                        for (int j = 0; j < imin(n, m); j++) {
                            for (i = j+1; i < m; i++) {
                                value[i] += fabsf(A[lda*j+i]);
                            }
                        }
                    }
                }
            }
        }
        break;
    }
}
