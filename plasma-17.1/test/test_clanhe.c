/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from test/test_zlanhe.c, normal z -> c, Sun Feb 20 23:48:08 2022
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
#include <math.h>

#include <omp.h>

#define COMPLEX

/***************************************************************************//**
 *
 * @brief Tests CLANHE.
 *
 * @param[in]  param - array of parameters
 * @param[out] info  - string of column labels or column values; length InfoLen
 *
 * If param is NULL and info is NULL,     print usage and return.
 * If param is NULL and info is non-NULL, set info to column labels and return.
 * If param is non-NULL and info is non-NULL, set info to column values
 * and run test.
 ******************************************************************************/
void test_clanhe(param_value_t param[], char *info)
{
    //================================================================
    // Print usage info or return column labels or values.
    //================================================================
    if (param == NULL) {
        if (info == NULL) {
            // Print usage info.
            print_usage(PARAM_NORM);
            print_usage(PARAM_UPLO);
            print_usage(PARAM_N);
            print_usage(PARAM_PADA);
            print_usage(PARAM_NB);
        }
        else {
            // Return column labels.
            snprintf(info, InfoLen,
                     "%*s %*s %*s %*s %*s",
                     InfoSpacing, "Norm",
                     InfoSpacing, "Uplo",
                     InfoSpacing, "N",
                     InfoSpacing, "PadA",
                     InfoSpacing, "NB");
        }
        return;
    }
    // Return column values.
    snprintf(info, InfoLen,
             "%*c %*c %*d %*d %*d",
             InfoSpacing, param[PARAM_NORM].c,
             InfoSpacing, param[PARAM_UPLO].c,
             InfoSpacing, param[PARAM_N].i,
             InfoSpacing, param[PARAM_PADA].i,
             InfoSpacing, param[PARAM_NB].i);

    //================================================================
    // Set parameters.
    //================================================================
    plasma_enum_t norm = plasma_norm_const(param[PARAM_NORM].c);
    plasma_enum_t uplo = plasma_uplo_const(param[PARAM_UPLO].c);

    int n = param[PARAM_N].i;

    int lda = imax(1, n + param[PARAM_PADA].i);

    int test = param[PARAM_TEST].c == 'y';
    float eps = LAPACKE_slamch('E');

    //================================================================
    // Set tuning parameters.
    //================================================================
    plasma_set(PlasmaNb, param[PARAM_NB].i);

    //================================================================
    // Allocate and initialize arrays.
    //================================================================
    plasma_complex32_t *A =
        (plasma_complex32_t*)malloc((size_t)lda*n*sizeof(plasma_complex32_t));
    assert(A != NULL);

    int seed[] = {0, 0, 0, 1};
    lapack_int retval;
    retval = LAPACKE_clarnv(1, seed, (size_t)lda*n, A);
    assert(retval == 0);

    plasma_complex32_t *Aref = NULL;
    float *work = NULL;
    if (test) {
        Aref = (plasma_complex32_t*)malloc(
            (size_t)lda*n*sizeof(plasma_complex32_t));
        assert(Aref != NULL);

        memcpy(Aref, A, (size_t)lda*n*sizeof(plasma_complex32_t));

        work = (float*)malloc(n*sizeof(float));
        assert(work != NULL);
    }

    //================================================================
    // Run and time PLASMA.
    //================================================================
    plasma_time_t start = omp_get_wtime();
    float value = plasma_clanhe(norm, uplo, n, A, lda);
    plasma_time_t stop = omp_get_wtime();
    plasma_time_t time = stop-start;

    param[PARAM_TIME].d = time;
    param[PARAM_GFLOPS].d = flops_clanhe(n, norm) / time / 1e9;

    //================================================================
    // Test results by comparing to a reference implementation.
    //================================================================
    if (test) {
        float valueRef =
            LAPACKE_clanhe_work(LAPACK_COL_MAJOR,
                                lapack_const(norm), lapack_const(uplo),
                                n, Aref, lda, work);

        // Calculate relative error
        float error = fabsf(value-valueRef);
        if (valueRef != 0)
            error /= valueRef;
        float tol = eps;
        float normalize = 1;
        switch (norm) {
            case PlasmaInfNorm:
                // Sum order on the line can differ
                normalize = n;
                break;

            case PlasmaOneNorm:
                // Sum order on the column can differ
                normalize = n;
                break;

            case PlasmaFrobeniusNorm:
                // Sum order on every element can differ
                normalize = n*n;
                break;
        }
        error /= normalize;
        param[PARAM_ERROR].d   = error;
        param[PARAM_SUCCESS].i = error < tol;
    }

    //================================================================
    // Free arrays.
    //================================================================
    free(A);
    if (test) {
        free(Aref);
        free(work);
    }
}
