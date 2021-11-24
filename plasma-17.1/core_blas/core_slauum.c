/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from core_blas/core_zlauum.c, normal z -> s, Mon Nov 22 19:22:13 2021
 *
 **/

#include "core_blas.h"
#include "plasma_types.h"
#include "core_lapack.h"

/***************************************************************************//**
 *
 * @ingroup core_lauum
 *
 *  Computes the product U * U^T or L^T * L, where the triangular
 *  factor U or L is stored in the upper or lower triangular part of
 *  the array A.
 *
 *  If uplo = 'U' or 'u' then the upper triangle of the result is stored,
 *  overwriting the factor U in A.
 *  If uplo = 'L' or 'l' then the lower triangle of the result is stored,
 *  overwriting the factor L in A.

 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = PlasmaUpper: Upper triangle of A is stored;
 *          = PlasmaLower: Lower triangle of A is stored.
 *
 *
 * @param[in] n
 *          The order of the matrix A. n >= 0.
 *
 * @param[in,out] A
 *          On entry, the triangular factor U or L.
 *          On exit, if uplo = 'U', the upper triangle of A is
 *          overwritten with the upper triangle of the product U * U^T;
 *          if uplo = 'L', the lower triangle of A is overwritten with
 *          the lower triangle of the product L^T * L.

 *
 * @param[in] lda
 *          The leading dimension of the array A. lda >= max(1,n).
 *
 * @param[out] info
 *          - 0 on successful exit
 *          - < 0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
int core_slauum(plasma_enum_t uplo,
                int n,
                float *A, int lda)
{
    return LAPACKE_slauum_work(LAPACK_COL_MAJOR,
                        lapack_const(uplo), n, A, lda);
}

/******************************************************************************/
void core_omp_slauum(plasma_enum_t uplo,
                     int n,
                     float *A, int lda,
                     plasma_sequence_t *sequence, plasma_request_t *request)
{
    #pragma omp task depend(inout:A[0:lda*n])
    {
        if (sequence->status == PlasmaSuccess) {
            int info = core_slauum(uplo, n, A, lda);
            if (info != PlasmaSuccess) {
                coreblas_error("core_slauum() failed");
                plasma_request_fail(sequence, request, PlasmaErrorInternal);
            }
        }
    }
}
