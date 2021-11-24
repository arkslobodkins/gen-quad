/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from test/test_zlag2c.c, mixed zc -> ds, Mon Nov 22 19:23:01 2021
 *
 **/

#include "core_blas.h"
#include "core_lapack.h"
#include "flops.h"
#include "plasma.h"
#include "test.h"

#include <assert.h>
#include <math.h>
#include <omp.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/***************************************************************************//**
 *
 * @brief Tests DLAG2S
 *
 * @param[in]  param - array of parameters
 * @param[out] info  - string of column labels or column values; length InfoLen
 *
 * If param is NULL and info is NULL,     print usage and return.
 * If param is NULL and info is non-NULL, set info to column labels and return.
 * If param is non-NULL and info is non-NULL, set info to column values
 * and run test.
 ******************************************************************************/
void test_dlag2s(param_value_t param[], char *info)
{
    //================================================================
    // Print usage info or return column labels or values
    //================================================================
    if (param == NULL) {
        if (info == NULL) {
            // Print usage info
            print_usage(PARAM_M);
            print_usage(PARAM_N);
            print_usage(PARAM_PADA);
            print_usage(PARAM_PADB);
            print_usage(PARAM_NB);
        }
        else {
            // Return column labels
            snprintf(info, InfoLen,
                     "%*s %*s %*s %*s %*s",
                     InfoSpacing, "m",
                     InfoSpacing, "n",
                     InfoSpacing, "PadA",
                     InfoSpacing, "PadB",
                     InfoSpacing, "nb");
        }
        return;
    }
    // Return column values
    snprintf(info, InfoLen,
             "%*d %*d %*d %*d %*d",
             InfoSpacing, param[PARAM_M].i,
             InfoSpacing, param[PARAM_N].i,
             InfoSpacing, param[PARAM_PADA].i,
             InfoSpacing, param[PARAM_PADB].i,
             InfoSpacing, param[PARAM_NB].i);

    //================================================================
    // Set parameters
    //================================================================
    int m = param[PARAM_M].i;
    int n = param[PARAM_N].i;

    int lda  = imax(1, m + param[PARAM_PADA].i);
    int ldas = imax(1, m + param[PARAM_PADB].i);

    int    test = param[PARAM_TEST].c == 'y';
    double eps  = LAPACKE_dlamch('E');

    //================================================================
    // Set tuning parameters
    //================================================================
    plasma_set(PlasmaNb, param[PARAM_NB].i);

    //================================================================
    // Allocate and initialize arrays
    //================================================================
    double *A =
        (double*)malloc((size_t)lda*n*sizeof(double));
    assert(A != NULL);

    float *As =
        (float*)malloc((size_t)ldas*n*sizeof(float));
    assert(As != NULL);

    int seed[] = {0, 0, 0, 1};
    lapack_int retval;
    retval = LAPACKE_dlarnv(1, seed, (size_t)lda*n, A);
    assert(retval == 0);

    float *AsRef = NULL;
    if (test) {
        AsRef = (float*)malloc(
            (size_t)ldas*n*sizeof(float));
        assert(AsRef != NULL);

        memcpy(AsRef, As, (size_t)ldas*n*sizeof(float));
    }

    //================================================================
    // Run and time PLASMA
    //================================================================
    plasma_time_t start = omp_get_wtime();

    retval = plasma_dlag2s(m, n, A, lda, As, ldas);

    plasma_time_t stop = omp_get_wtime();

    param[PARAM_GFLOPS].d  = 0.0;

    if (retval != PlasmaSuccess) {
        plasma_error("plasma_dlag2s() failed");
        param[PARAM_TIME].d    = 0.0;
        param[PARAM_ERROR].d   = 1.0;
        param[PARAM_SUCCESS].i = false;
        return;
    }
    else {
        param[PARAM_TIME].d = stop-start;
    }

    //================================================================
    // Test results by comparing to result of LAPACK function
    //================================================================
    if (test) {
        // Calculate relative error |As_ref - As|_F / |As_ref|_F < 3*eps
        // Using 3*eps covers complex arithmetic

        lapack_int mtrxLayout = LAPACK_COL_MAJOR;

        retval = LAPACKE_dlag2s_work(mtrxLayout, m, n, A, lda, AsRef, ldas);

        if (retval != PlasmaSuccess) {
            coreblas_error("LAPACKE_dlag2s_work() failed");
            param[PARAM_ERROR].d   = 1.0;
            param[PARAM_SUCCESS].i = false;
            return;
        }

        float work[1];

        // Calculate Frobenius norm of reference result As_ref
        double AsNormRef = LAPACKE_slange_work(mtrxLayout, 'F',
                                               ldas, n, AsRef, ldas, work);

        // Calculate difference As_ref-As
        float cmone = -1.0;
        cblas_saxpy((size_t)ldas*n, (cmone), As, 1, AsRef, 1);

        // Calculate Frobenius norm of As_ref-As
        double AsNormDiff = LAPACKE_slange_work(mtrxLayout, 'F',
                                                ldas, n, AsRef, ldas, work);

        // Calculate relative error |As_ref-As|_F / |As_ref|_F
        double error = AsNormDiff/AsNormRef;

        param[PARAM_ERROR].d   = error;
        param[PARAM_SUCCESS].i = error < 3*eps;
    }

    //================================================================
    // Free arrays
    //================================================================
    free(A);
    free(As);

    if (test)
        free(AsRef);
}
