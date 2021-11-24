/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from test/test_zsyrk.c, normal z -> s, Mon Nov 22 19:22:54 2021
 *
 **/
#include "test.h"
#include "flops.h"
#include "core_lapack.h"
#include "plasma.h"

#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <omp.h>

#define REAL

/***************************************************************************//**
 *
 * @brief Tests SSYRK.
 *
 * @param[in]  param - array of parameters
 * @param[out] info  - string of column labels or column values; length InfoLen
 *
 * If param is NULL and info is NULL,     print usage and return.
 * If param is NULL and info is non-NULL, set info to column labels and return.
 * If param is non-NULL and info is non-NULL, set info to column values
 * and run test.
 ******************************************************************************/
void test_ssyrk(param_value_t param[], char *info)
{
    //================================================================
    // Print usage info or return column labels or values.
    //================================================================
    if (param == NULL) {
        if (info == NULL) {
            // Print usage info.
            print_usage(PARAM_UPLO);
            print_usage(PARAM_TRANS);
            print_usage(PARAM_N);
            print_usage(PARAM_K);
            print_usage(PARAM_ALPHA);
            print_usage(PARAM_BETA);
            print_usage(PARAM_PADA);
            print_usage(PARAM_PADC);
            print_usage(PARAM_NB);
        }
        else {
            // Return column labels.
            snprintf(info, InfoLen,
                     "%*s %*s %*s %*s %*s %*s %*s %*s %*s",
                     InfoSpacing, "Uplo",
                     InfoSpacing, "Trans",
                     InfoSpacing, "N",
                     InfoSpacing, "K",
                     InfoSpacing, "alpha",
                     InfoSpacing, "beta",
                     InfoSpacing, "PadA",
                     InfoSpacing, "PadC",
                     InfoSpacing, "NB");
        }
        return;
    }
    // Return column values.
    snprintf(info, InfoLen,
             "%*c %*c %*d %*d %*.4f %*.4f %*d %*d %*d",
             InfoSpacing, param[PARAM_UPLO].c,
             InfoSpacing, param[PARAM_TRANS].c,
             InfoSpacing, param[PARAM_N].i,
             InfoSpacing, param[PARAM_K].i,
             InfoSpacing, creal(param[PARAM_ALPHA].z),
             InfoSpacing, creal(param[PARAM_BETA].z),
             InfoSpacing, param[PARAM_PADA].i,
             InfoSpacing, param[PARAM_PADC].i,
             InfoSpacing, param[PARAM_NB].i);

    //================================================================
    // Set parameters.
    //================================================================
    plasma_enum_t uplo = plasma_uplo_const(param[PARAM_UPLO].c);
    plasma_enum_t trans = plasma_trans_const(param[PARAM_TRANS].c);

    int n = param[PARAM_N].i;
    int k = param[PARAM_K].i;

    int Am, An;
    int Cm, Cn;

    if (trans == PlasmaNoTrans) {
        Am = n;
        An = k;
    }
    else {
        Am = k;
        An = n;
    }
    Cm = n;
    Cn = n;

    int lda = imax(1, Am + param[PARAM_PADA].i);
    int ldc = imax(1, Cm + param[PARAM_PADC].i);

    int test = param[PARAM_TEST].c == 'y';
    float eps = LAPACKE_slamch('E');

#ifdef COMPLEX
    float alpha = param[PARAM_ALPHA].z;
    float beta  = param[PARAM_BETA].z;
#else
    float alpha = creal(param[PARAM_ALPHA].z);
    float beta  = creal(param[PARAM_BETA].z);
#endif

    //================================================================
    // Set tuning parameters.
    //================================================================
    plasma_set(PlasmaNb, param[PARAM_NB].i);

    //================================================================
    // Allocate and initialize arrays.
    //================================================================
    float *A =
        (float*)malloc((size_t)lda*An*sizeof(float));
    assert(A != NULL);

    float *C =
        (float*)malloc((size_t)ldc*Cn*sizeof(float));
    assert(C != NULL);

    int seed[] = {0, 0, 0, 1};
    lapack_int retval;
    retval = LAPACKE_slarnv(1, seed, (size_t)lda*An, A);
    assert(retval == 0);

    retval = LAPACKE_slarnv(1, seed, (size_t)ldc*Cn, C);
    assert(retval == 0);

    float *Cref = NULL;
    if (test) {
        Cref = (float*)malloc(
            (size_t)ldc*Cn*sizeof(float));
        assert(Cref != NULL);

        memcpy(Cref, C, (size_t)ldc*Cn*sizeof(float));
    }

    //================================================================
    // Run and time PLASMA.
    //================================================================
    plasma_time_t start = omp_get_wtime();

    plasma_ssyrk(
        uplo, trans,
        n, k,
        alpha, A, lda,
        beta, C, ldc);

    plasma_time_t stop = omp_get_wtime();
    plasma_time_t time = stop-start;

    param[PARAM_TIME].d = time;
    param[PARAM_GFLOPS].d = flops_ssyrk(n, k) / time / 1e9;

    //================================================================
    // Test results by comparing to a reference implementation.
    //================================================================
    if (test) {
        // see comments in test_sgemm.c
        char uplo_ = param[PARAM_UPLO].c;
        float work[1];
        float Anorm = LAPACKE_slange_work(
                           LAPACK_COL_MAJOR, 'F', Am, An, A, lda, work);
        float Cnorm = LAPACKE_slansy_work(
                           LAPACK_COL_MAJOR, 'F', uplo_, Cn, Cref, ldc, work);

        cblas_ssyrk(
            CblasColMajor,
            (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
            n, k,
            (alpha), A, lda,
            (beta), Cref, ldc);

        float zmone = -1.0;
        cblas_saxpy((size_t)ldc*Cn, (zmone), Cref, 1, C, 1);

        float error = LAPACKE_slansy_work(
                           LAPACK_COL_MAJOR, 'F', uplo_, Cn, C,    ldc, work);
        float normalize = sqrtf((float)k+2) * fabsf(alpha) * Anorm * Anorm
                         + 2 * fabsf(beta) * Cnorm;
        if (normalize != 0)
            error /= normalize;

        param[PARAM_ERROR].d = error;
        param[PARAM_SUCCESS].i = error < 3*eps;
    }

    //================================================================
    // Free arrays.
    //================================================================
    free(A);
    free(C);
    if (test)
        free(Cref);
}
