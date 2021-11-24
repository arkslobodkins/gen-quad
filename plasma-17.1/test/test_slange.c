/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from test/test_zlange.c, normal z -> s, Mon Nov 22 19:22:52 2021
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

#define REAL

/***************************************************************************//**
 *
 * @brief Tests SLANGE.
 *
 * @param[in]  param - array of parameters
 * @param[out] info  - string of column labels or column values; length InfoLen
 *
 * If param is NULL and info is NULL,     print usage and return.
 * If param is NULL and info is non-NULL, set info to column labels and return.
 * If param is non-NULL and info is non-NULL, set info to column values
 * and run test.
 ******************************************************************************/
void test_slange(param_value_t param[], char *info)
{
    //================================================================
    // Print usage info or return column labels or values.
    //================================================================
    if (param == NULL) {
        if (info == NULL) {
            // Print usage info.
            print_usage(PARAM_NORM);
            print_usage(PARAM_M);
            print_usage(PARAM_N);
            print_usage(PARAM_PADA);
            print_usage(PARAM_NB);
        }
        else {
            // Return column labels.
            snprintf(info, InfoLen,
                     "%*s %*s %*s %*s %*s",
                     InfoSpacing, "Norm",
                     InfoSpacing, "M",
                     InfoSpacing, "N",
                     InfoSpacing, "PadA",
                     InfoSpacing, "NB");
        }
        return;
    }
    // Return column values.
    snprintf(info, InfoLen,
             "%*c %*d %*d %*d %*d",
             InfoSpacing, param[PARAM_NORM].c,
             InfoSpacing, param[PARAM_M].i,
             InfoSpacing, param[PARAM_N].i,
             InfoSpacing, param[PARAM_PADA].i,
             InfoSpacing, param[PARAM_NB].i);

    //================================================================
    // Set parameters.
    //================================================================
    plasma_enum_t norm = plasma_norm_const(param[PARAM_NORM].c);

    int m = param[PARAM_M].i;
    int n = param[PARAM_N].i;

    int lda = imax(1, m + param[PARAM_PADA].i);

    int test = param[PARAM_TEST].c == 'y';
    float eps = LAPACKE_slamch('E');

    //================================================================
    // Set tuning parameters.
    //================================================================
    plasma_set(PlasmaNb, param[PARAM_NB].i);

    //================================================================
    // Allocate and initialize arrays.
    //================================================================
    float *A =
        (float*)malloc((size_t)lda*n*sizeof(float));
    assert(A != NULL);

    int seed[] = {0, 0, 0, 1};
    lapack_int retval;
    retval = LAPACKE_slarnv(1, seed, (size_t)lda*n, A);
    assert(retval == 0);

    float *Aref = NULL;
    if (test) {
        Aref = (float*)malloc(
            (size_t)lda*n*sizeof(float));
        assert(Aref != NULL);

        memcpy(Aref, A, (size_t)lda*n*sizeof(float));
    }

    //================================================================
    // Run and time PLASMA.
    //================================================================
    plasma_time_t start = omp_get_wtime();
    float value = plasma_slange(norm, m, n, A, lda);
    plasma_time_t stop = omp_get_wtime();
    plasma_time_t time = stop-start;

    param[PARAM_TIME].d = time;
    param[PARAM_GFLOPS].d = flops_slange(m, n, norm) / time / 1e9;

    //================================================================
    // Test results by comparing to a reference implementation.
    //================================================================
    if (test) {
        float valueRef =
            LAPACKE_slange(LAPACK_COL_MAJOR, lapack_const(norm),
                           m, n, Aref, lda);

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
                normalize = m;
                break;

            case PlasmaFrobeniusNorm:
                // Sum order on every element can differ
                normalize = m*n;
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
    if (test)
        free(Aref);
}
