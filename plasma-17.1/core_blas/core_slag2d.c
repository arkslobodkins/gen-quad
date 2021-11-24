/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from core_blas/core_clag2z.c, mixed zc -> ds, Mon Nov 22 19:22:22 2021
 *
 **/

#include "core_blas.h"
#include "core_lapack.h"
#include "plasma_types.h"

/***************************************************************************//**
 *
 * @ingroup core_lag2
 *
 *  Converts m-by-n matrix A from single complex to double complex precision.
 *
 *******************************************************************************
 *
 * @param[in] m
 *          The number of rows of the matrix As.
 *          m >= 0.
 *
 * @param[in] n
 *          The number of columns of the matrix As.
 *          n >= 0.
 *
 * @param[in] As
 *          The ldas-by-n matrix in single complex precision to convert.
 *
 * @param[in] ldas
 *          The leading dimension of the matrix As.
 *          ldas >= max(1,m).
 *
 * @param[out] A
 *          On exit, the converted lda-by-n matrix in double complex precision.
 *
 * @param[in] lda
 *          The leading dimension of the matrix A.
 *          lda >= max(1,m).
 *
 ******************************************************************************/
void core_slag2d(int m, int n,
                 float *As, int ldas,
                 double *A,  int lda)
{
    LAPACKE_slag2d_work(LAPACK_COL_MAJOR, m, n, As, ldas, A, lda);
}

/******************************************************************************/
void core_omp_slag2d(int m, int n,
                     float *As, int ldas,
                     double *A,  int lda,
                     plasma_sequence_t *sequence, plasma_request_t *request)
{
    #pragma omp task depend(in:As[0:ldas*n]) \
                     depend(out:A[0:lda*n])
    {
        if (sequence->status == PlasmaSuccess)
            core_slag2d(m, n, As, ldas, A, lda);
    }
}
