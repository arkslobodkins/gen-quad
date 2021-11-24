/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from test/test_zlaswp.c, normal z -> s, Mon Nov 22 19:23:13 2021
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
 * @brief Tests SLASWP
 *
 * @param[in]  param - array of parameters
 * @param[out] info  - string of column labels or column values; length InfoLen
 *
 * If param is NULL and info is NULL,     print usage and return.
 * If param is NULL and info is non-NULL, set info to column labels and return.
 * If param is non-NULL and info is non-NULL, set info to column values
 * and run test.
 ******************************************************************************/
void test_slaswp(param_value_t param[], char *info)
{
    //================================================================
    // Print usage info or return column labels or values
    //================================================================
    if (param == NULL) {
        if (info == NULL) {
            // Print usage info
            print_usage(PARAM_COLROW);
            print_usage(PARAM_M);
            print_usage(PARAM_N);
            print_usage(PARAM_PADA);
            print_usage(PARAM_NB);
            print_usage(PARAM_INCX);
        }
        else {
            // Return column labels
            snprintf(info, InfoLen,
                     "%*s %*s %*s %*s %*s %*s",
                     InfoSpacing, "ColRow",
                     InfoSpacing, "m",
                     InfoSpacing, "n",
                     InfoSpacing, "PadA",
                     InfoSpacing, "nb",
                     InfoSpacing, "incx");
        }
        return;
    }
    // Return column values
    snprintf(info, InfoLen,
             "%*c %*d %*d %*d %*d %*d",
             InfoSpacing, param[PARAM_COLROW].c,
             InfoSpacing, param[PARAM_M].i,
             InfoSpacing, param[PARAM_N].i,
             InfoSpacing, param[PARAM_PADA].i,
             InfoSpacing, param[PARAM_NB].i,
             InfoSpacing, param[PARAM_INCX].i);

    //================================================================
    // Set parameters
    //================================================================
    plasma_enum_t colrow = plasma_storev_const(param[PARAM_COLROW].c);
    int incx = param[PARAM_INCX].i;

    int m = param[PARAM_M].i;
    int n = param[PARAM_N].i;

    int lda = imax(1, m + param[PARAM_PADA].i);

    int test = param[PARAM_TEST].c == 'y';
    float tol = param[PARAM_TOL].d * LAPACKE_slamch('E');

    //================================================================
    // Set tuning parameters
    //================================================================
    plasma_set(PlasmaNb, param[PARAM_NB].i);

    //================================================================
    // Allocate and initialize arrays
    //================================================================
    float *A =
        (float*)malloc((size_t)lda*n*sizeof(float));
    assert(A != NULL);

    int *ipiv;
    if (colrow == PlasmaRowwise) {
        ipiv = (int*)malloc((size_t)m*sizeof(int));

        assert(ipiv != NULL);
        for (int i = 0; i < m; i++)
            ipiv[i] = rand()%m + 1;
    }
    else {
        ipiv = (int*)malloc((size_t)n*sizeof(int));
        assert(ipiv != NULL);

        for (int j = 0; j < n; j++)
            ipiv[j] = rand()%n + 1;
    }

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
    // Run and time PLASMA
    //================================================================
    plasma_time_t start = omp_get_wtime();
    retval = plasma_slaswp(colrow, m, n, A, lda, ipiv, incx);
    plasma_time_t stop = omp_get_wtime();

    param[PARAM_TIME].d = stop-start;
    param[PARAM_GFLOPS].d = 0.0;

    //================================================================
    // Test results by comparing to result of core_slacpy function
    //================================================================
    if (test) {
        if (colrow == PlasmaRowwise)
            LAPACKE_slaswp(LAPACK_COL_MAJOR,
                           n,
                           Aref, lda,
                           1, m, ipiv, incx);
        else
            LAPACKE_slaswp(LAPACK_ROW_MAJOR,
                           m,
                           Aref, lda,
                           1, n, ipiv, incx);

        float zmone = -1.0;
        cblas_saxpy((size_t)lda*n, (zmone), Aref, 1, A, 1);

        float work[1];
        float Anorm = LAPACKE_slange_work(
            LAPACK_COL_MAJOR, 'F', m, n, Aref, lda, work);

        float error = LAPACKE_slange_work(
            LAPACK_COL_MAJOR, 'F', m, n, A, lda, work);

        if (Anorm != 0)
            error /= Anorm;

        param[PARAM_ERROR].d = error;
        param[PARAM_SUCCESS].i = error < tol;
    }

    //================================================================
    // Free arrays
    //================================================================
    free(A);
    if (test)
        free(Aref);
}
