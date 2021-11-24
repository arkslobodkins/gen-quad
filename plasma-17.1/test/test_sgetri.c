/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from test/test_zgetri.c, normal z -> s, Mon Nov 22 19:22:46 2021
 *
 **/
#include "test.h"
#include "flops.h"
#include "core_blas.h"
#include "core_lapack.h"
#include "plasma.h"

#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <omp.h>

#define REAL

#define A(i_, j_) A[(i_) + (size_t)lda*(j_)]

/***************************************************************************//**
 *
 * @brief Tests strtri.
 *
 * @param[in]  param - array of parameters
 * @param[out] info  - string of column labels or column values; length InfoLen
 *
 * If param is NULL and info is NULL,     print usage and return.
 * If param is NULL and info is non-NULL, set info to column labels and return.
 * If param is non-NULL and info is non-NULL, set info to column values
 * and run test.
 ******************************************************************************/
void test_sgetri(param_value_t param[], char *info)
{
    //================================================================
    // Print usage info or return column labels or values.
    //================================================================
    if (param == NULL) {
        if (info == NULL) {
            // Print usage info.
            print_usage(PARAM_N);
            print_usage(PARAM_PADA);
            print_usage(PARAM_NB);
            print_usage(PARAM_ZEROCOL);
        }
        else {
            // Return column labels.
            snprintf(info, InfoLen,
                     "%*s %*s %*s %*s",
                     InfoSpacing, "N",
                     InfoSpacing, "PadA",
                     InfoSpacing, "NB",
                     InfoSpacing, "ZeroCol");
        }
        return;
    }
    // Return column values.
    snprintf(info, InfoLen,
             "%*d %*d %*d %*d",
             InfoSpacing, param[PARAM_N].i,
             InfoSpacing, param[PARAM_PADA].i,
             InfoSpacing, param[PARAM_NB].i,
             InfoSpacing, param[PARAM_ZEROCOL].i);

    //================================================================
    // Set parameters.
    //================================================================
    int n = param[PARAM_N].i;
    int lda = imax(1, n + param[PARAM_PADA].i);

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
        (float*)malloc((size_t)lda*n*sizeof(float));
    assert(A != NULL);

    int *ipiv;
    ipiv = (int*)malloc((size_t)n*sizeof(int));
    assert(ipiv != NULL);

    int seed[] = {0, 0, 0, 1};
    lapack_int retval;
    retval = LAPACKE_slarnv(1, seed, (size_t)lda*n, A);
    assert(retval == 0);

    // Take LU decomposition
    LAPACKE_sgetrf(LAPACK_COL_MAJOR, n, n, A, lda, ipiv);

    int zerocol = param[PARAM_ZEROCOL].i;
    if (zerocol >= 0 && zerocol < n)
        memset(&A[zerocol*lda], 0, n*sizeof(float));

    float *Aref;
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
    int plainfo = plasma_sgetri(n, A, lda, ipiv);
    plasma_time_t stop = omp_get_wtime();
    plasma_time_t time = stop-start;

    param[PARAM_TIME].d = time;
    param[PARAM_GFLOPS].d = flops_strtri(n) / time / 1e9;

    //================================================================
    // Test results by checking the residual
    // ||B - A|| / (||A||)
    //================================================================
    if (test) {
        float zmone = -1.0;

        // norm(A)
        float temp;
        float Anorm = LAPACKE_slange_work(
               LAPACK_COL_MAJOR, 'F', n, n, Aref, lda, &temp);

        int lwork = n;
        float *work =
            (float*)malloc(lwork * sizeof(float));

        // B = inv(A)
        int lapinfo = LAPACKE_sgetri_work(CblasColMajor,
                                          n, Aref, lda,
                                          ipiv, work, lwork);
        if (lapinfo == 0) {
            // norm(A^{-1})
            float Inorm = LAPACKE_slange_work(
                   LAPACK_COL_MAJOR, 'F', n, n, Aref, lda, &temp);

            // A -= Aref
            cblas_saxpy((size_t)lda*n, (zmone), Aref, 1, A, 1);

            float error = LAPACKE_slange_work(
                               LAPACK_COL_MAJOR, 'F', lda, n, A, lda, &temp);
            if (Anorm*Inorm != 0)
                error /= (Anorm * Inorm);

            param[PARAM_ERROR].d = error;
            param[PARAM_SUCCESS].i = error < tol;
        }
        else {
            if (plainfo == lapinfo) {
                param[PARAM_ERROR].d = 0.0;
                param[PARAM_SUCCESS].i = 1;
            }
            else {
                param[PARAM_ERROR].d = INFINITY;
                param[PARAM_SUCCESS].i = 0;
            }
        }
        free(work);
    }

    //================================================================
    // Free arrays.
    //================================================================
    free(A);
    free(ipiv);
    if (test)
        free(Aref);
}
