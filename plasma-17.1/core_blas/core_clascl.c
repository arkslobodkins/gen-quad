/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from core_blas/core_zlascl.c, normal z -> c, Tue Feb  8 19:16:20 2022
 *
 **/

#include "core_blas.h"
#include "plasma_types.h"
#include "core_lapack.h"

/******************************************************************************/
void core_clascl(plasma_enum_t uplo,
                 float cfrom, float cto,
                 int m, int n,
                 plasma_complex32_t *A, int lda)
{
    // LAPACKE_clascl is not available in LAPACKE < 3.6.0
    int kl;
    int ku;
    int info;
    char type = lapack_const(uplo);
    LAPACK_clascl(&type,
                  &kl, &ku,
                  &cfrom, &cto,
                  &m, &n,
                  A, &lda, &info);
}

/******************************************************************************/
void core_omp_clascl(plasma_enum_t uplo,
                     float cfrom, float cto,
                     int m, int n,
                     plasma_complex32_t *A, int lda,
                     plasma_sequence_t *sequence, plasma_request_t *request)
{
    #pragma omp task depend(inout:A[0:lda*n])
    {
        if (sequence->status == PlasmaSuccess)
            core_clascl(uplo,
                        cfrom, cto,
                        m, n,
                        A, lda);
    }
}
