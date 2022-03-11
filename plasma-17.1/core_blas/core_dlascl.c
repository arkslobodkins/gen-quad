/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from core_blas/core_zlascl.c, normal z -> d, Thu Mar 10 18:58:33 2022
 *
 **/

#include "core_blas.h"
#include "plasma_types.h"
#include "core_lapack.h"

/******************************************************************************/
void core_dlascl(plasma_enum_t uplo,
                 double cfrom, double cto,
                 int m, int n,
                 double *A, int lda)
{
    // LAPACKE_dlascl is not available in LAPACKE < 3.6.0
    int kl;
    int ku;
    int info;
    char type = lapack_const(uplo);
    LAPACK_dlascl(&type,
                  &kl, &ku,
                  &cfrom, &cto,
                  &m, &n,
                  A, &lda, &info);
}

/******************************************************************************/
void core_omp_dlascl(plasma_enum_t uplo,
                     double cfrom, double cto,
                     int m, int n,
                     double *A, int lda,
                     plasma_sequence_t *sequence, plasma_request_t *request)
{
    #pragma omp task depend(inout:A[0:lda*n])
    {
        if (sequence->status == PlasmaSuccess)
            core_dlascl(uplo,
                        cfrom, cto,
                        m, n,
                        A, lda);
    }
}
