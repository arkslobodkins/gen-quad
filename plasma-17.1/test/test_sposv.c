/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from test/test_zposv.c, normal z -> s, Mon Nov 22 19:22:55 2021
 *
 **/
#include "test.h"
#include "flops.h"
#include "core_blas.h"
#include "core_lapack.h"
#include "plasma.h"

#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <omp.h>

#define REAL

#define A(i_, j_) A[(i_) + (size_t)lda*(j_)]

/***************************************************************************//**
 *
 * @brief Tests SPOSV.
 *
 * @param[in]  param - array of parameters
 * @param[out] info  - string of column labels or column values; length InfoLen
 *
 * If param is NULL and info is NULL,     print usage and return.
 * If param is NULL and info is non-NULL, set info to column labels and return.
 * If param is non-NULL and info is non-NULL, set info to column values
 * and run test.
 ******************************************************************************/
void test_sposv(param_value_t param[], char *info)
{
    //================================================================
    // Print usage info or return column labels or values.
    //================================================================
    if (param == NULL) {
        if (info == NULL) {
            // Print usage info.
            print_usage(PARAM_UPLO);
            print_usage(PARAM_N);
            print_usage(PARAM_NRHS);
            print_usage(PARAM_PADA);
            print_usage(PARAM_PADB);
            print_usage(PARAM_NB);
        }
        else {
            // Return column labels.
            snprintf(info, InfoLen,
                     "%*s %*s %*s %*s %*s %*s",
                     InfoSpacing, "Uplo",
                     InfoSpacing, "N",
                     InfoSpacing, "NRHS",
                     InfoSpacing, "PadA",
                     InfoSpacing, "PadB",
                     InfoSpacing, "NB");
        }
        return;
    }
    // Return column values.
    snprintf(info, InfoLen,
             "%*c %*d %*d %*d %*d %*d",
             InfoSpacing, param[PARAM_UPLO].c,
             InfoSpacing, param[PARAM_N].i,
             InfoSpacing, param[PARAM_NRHS].i,
             InfoSpacing, param[PARAM_PADA].i,
             InfoSpacing, param[PARAM_PADB].i,
             InfoSpacing, param[PARAM_NB].i);

    //================================================================
    // Set parameters.
    //================================================================
    plasma_enum_t uplo = plasma_uplo_const(param[PARAM_UPLO].c);

    int n = param[PARAM_N].i;
    int nrhs = param[PARAM_NRHS].i;

    int lda = imax(1, n + param[PARAM_PADA].i);
    int ldb = imax(1, n + param[PARAM_PADB].i);

    int test = param[PARAM_TEST].c == 'y';
    float tol = param[PARAM_TOL].d * LAPACKE_slamch('E');

    //================================================================
    // Set tuning parameters.
    //================================================================
    plasma_set(PlasmaNb, param[PARAM_NB].i);

    //================================================================
    // Allocate and initialize arrays.
    //================================================================
    float *A =
        (float*)malloc((size_t)lda*n
                                    *sizeof(float));
    assert(A != NULL);

    float *B =
        (float*)malloc((size_t)ldb*nrhs
                                    *sizeof(float));
    assert(B != NULL);

    int seed[] = {0, 0, 0, 1};
    lapack_int retval;
    retval = LAPACKE_slarnv(1, seed, (size_t)lda*n, A);
    assert(retval == 0);

    retval = LAPACKE_slarnv(1, seed, (size_t)ldb*nrhs, B);
    assert(retval == 0);

    //================================================================
    // Make the A matrix symmetric/symmetric positive definite.
    // It increases diagonal by n, and makes it real.
    // It sets Aji = ( Aij ) for j < i, that is, copy lower
    // triangle to upper triangle.
    //================================================================
    for (int i = 0; i < n; ++i) {
        A(i,i) = creal(A(i,i)) + n;
        for (int j = 0; j < i; ++j) {
            A(j,i) = (A(i,j));
        }
    }

    float *Aref = NULL;
    float *Bref = NULL;
    float *work = NULL;
    if (test) {
        Aref = (float*)malloc(
            (size_t)lda*n*sizeof(float));
        assert(Aref != NULL);

        Bref = (float*)malloc(
            (size_t)ldb*nrhs*sizeof(float));
        assert(Bref != NULL);

        memcpy(Aref, A, (size_t)lda*n*sizeof(float));
        memcpy(Bref, B, (size_t)ldb*nrhs*sizeof(float));
    }

    //================================================================
    // Run and time PLASMA.
    //================================================================
    plasma_time_t start = omp_get_wtime();
    plasma_sposv(uplo, n, nrhs, A, lda, B, ldb);
    plasma_time_t stop = omp_get_wtime();
    plasma_time_t time = stop-start;

    float flops = flops_spotrf(n) + flops_spotrs(n, nrhs);
    param[PARAM_TIME].d = time;
    param[PARAM_GFLOPS].d = flops / time / 1e9;

    //================================================================
    // Test results by checking the residual
    //
    //                      || B - AX ||_I
    //                --------------------------- < epsilon
    //                 || A ||_I * || X ||_I * N
    //
    //================================================================
    if (test) {
        float zone  =  1.0;
        float zmone = -1.0;

        work = (float*)malloc((size_t)n*sizeof(float));
        assert(work != NULL);

        float Anorm = LAPACKE_slansy_work(
            LAPACK_COL_MAJOR, 'I', lapack_const(uplo), n, Aref, lda, work);
        float Xnorm = LAPACKE_slange_work(
            LAPACK_COL_MAJOR, 'I', n, nrhs, B, ldb, work);

        // Bref -= Aref*B
        cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, nrhs, n,
                    (zmone), Aref, lda,
                                        B,    ldb,
                    (zone),  Bref, ldb);

        float Rnorm = LAPACKE_slange_work(
            LAPACK_COL_MAJOR, 'I', n, nrhs, Bref, ldb, work);
        float residual = Rnorm/(n*Anorm*Xnorm);

        param[PARAM_ERROR].d = residual;
        param[PARAM_SUCCESS].i = residual < tol;
    }

    //================================================================
    // Free arrays.
    //================================================================
    free(A);
    free(B);
    if (test) {
        free(Aref);
        free(Bref);
        free(work);
    }
}
