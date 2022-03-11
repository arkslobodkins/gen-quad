/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from core_blas/core_zlacpy.c, normal z -> c, Thu Mar 10 18:58:38 2022
 *
 **/

#include "core_blas.h"
#include "plasma_types.h"
#include "core_lapack.h"

/***************************************************************************//**
 *
 * @ingroup core_lacpy
 *
 *  Copies all or part of a two-dimensional matrix A to another matrix B.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          - PlasmaGeneral: entire A,
 *          - PlasmaUpper:   upper triangle,
 *          - PlasmaLower:   lower triangle.
 *
 * @param[in] m
 *          The number of rows of the matrices A and B.
 *          m >= 0.
 *
 * @param[in] n
 *          The number of columns of the matrices A and B.
 *          n >= 0.
 *
 * @param[in] A
 *          The m-by-n matrix to copy.
 *
 * @param[in] lda
 *          The leading dimension of the array A.
 *          lda >= max(1,m).
 *
 * @param[out] B
 *          The m-by-n copy of the matrix A.
 *          On exit, B = A ONLY in the locations specified by uplo.
 *
 * @param[in] ldb
 *          The leading dimension of the array B.
 *          ldb >= max(1,m).
 *
 ******************************************************************************/
void core_clacpy(plasma_enum_t uplo,
                 int m, int n,
                 const plasma_complex32_t *A, int lda,
                       plasma_complex32_t *B, int ldb)
{
    LAPACKE_clacpy_work(LAPACK_COL_MAJOR,
                        lapack_const(uplo),
                        m, n,
                        A, lda,
                        B, ldb);
}

/******************************************************************************/
void core_omp_clacpy(plasma_enum_t uplo,
                     int m, int n,
                     const plasma_complex32_t *A, int lda,
                           plasma_complex32_t *B, int ldb,
                     plasma_sequence_t *sequence, plasma_request_t *request)
{
    #pragma omp task depend(in:A[0:lda*n]) \
                     depend(out:B[0:ldb*n])
    {
        if (sequence->status == PlasmaSuccess)
            core_clacpy(uplo,
                        m, n,
                        A, lda,
                        B, ldb);
    }
}
