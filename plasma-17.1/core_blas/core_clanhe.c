/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from core_blas/core_zlanhe.c, normal z -> c, Mon Nov 22 19:22:16 2021
 *
 **/

#include "core_blas.h"
#include "plasma_types.h"
#include "core_lapack.h"

#include <math.h>

/******************************************************************************/
void core_clanhe(plasma_enum_t norm, plasma_enum_t uplo,
                 int n,
                 const plasma_complex32_t *A, int lda,
                 float *work, float *value)
{
    *value = LAPACKE_clanhe_work(LAPACK_COL_MAJOR,
                                 lapack_const(norm),
                                 lapack_const(uplo),
                                 n, A, lda, work);
}

/******************************************************************************/
void core_omp_clanhe(plasma_enum_t norm, plasma_enum_t uplo,
                     int n,
                     const plasma_complex32_t *A, int lda,
                     float *work, float *value,
                     plasma_sequence_t *sequence, plasma_request_t *request)
{
    #pragma omp task depend(in:A[0:lda*n]) \
                     depend(out:value[0:1])
    {
        if (sequence->status == PlasmaSuccess)
            core_clanhe(norm, uplo, n, A, lda, work, value);
    }
}

/******************************************************************************/
void core_omp_clanhe_aux(plasma_enum_t norm, plasma_enum_t uplo,
                         int n,
                         const plasma_complex32_t *A, int lda,
                         float *value,
                         plasma_sequence_t *sequence, plasma_request_t *request)
{
    switch (norm) {
    case PlasmaOneNorm:
    case PlasmaInfNorm:
        #pragma omp task depend(in:A[0:lda*n]) \
                         depend(out:value[0:n])
        {
            if (sequence->status == PlasmaSuccess) {
                if (uplo == PlasmaUpper) {
                    for (int i = 0; i < n; i++)
                        value[i] = 0.0;

                    for (int j = 0; j < n; j++) {
                        for (int i = 0; i < j; i++) {
                            value[i] += cabsf(A[lda*j+i]);
                            value[j] += cabsf(A[lda*j+i]);
                        }
                        value[j] += fabsf(creal(A[lda*j+j]));
                    }
                }
                else { // PlasmaLower
                    for (int i = 0; i < n; i++)
                        value[i] = 0.0;

                    for (int j = 0; j < n; j++) {
                        value[j] += fabsf(creal(A[lda*j+j]));
                        for (int i = j+1; i < n; i++) {
                            value[i] += cabsf(A[lda*j+i]);
                            value[j] += cabsf(A[lda*j+i]);
                        }
                    }
                }
            }
        }
        break;
    }
}
