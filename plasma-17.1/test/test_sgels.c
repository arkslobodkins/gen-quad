/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from test/test_zgels.c, normal z -> s, Mon Nov 22 19:22:57 2021
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

#include <omp.h>

#define REAL

/***************************************************************************//**
 *
 * @brief Tests SGELS.
 *
 * @param[in]  param - array of parameters
 * @param[out] info  - string of column labels or column values; length InfoLen
 *
 * If param is NULL and info is NULL,     print usage and return.
 * If param is NULL and info is non-NULL, set info to column labels and return.
 * If param is non-NULL and info is non-NULL, set info to column values
 * and run test.
 ******************************************************************************/
void test_sgels(param_value_t param[], char *info)
{
    //================================================================
    // Print usage info or return column labels or values.
    //================================================================
    if (param == NULL) {
        if (info == NULL) {
            // Print usage info.
            print_usage(PARAM_M);
            print_usage(PARAM_N);
            print_usage(PARAM_NRHS);
            print_usage(PARAM_PADA);
            print_usage(PARAM_PADB);
            print_usage(PARAM_NB);
            print_usage(PARAM_IB);
            print_usage(PARAM_HMODE);
        }
        else {
            // Return column labels.
            snprintf(info, InfoLen,
                     "%*s %*s %*s %*s %*s %*s %*s %*s",
                     InfoSpacing, "M",
                     InfoSpacing, "N",
                     InfoSpacing, "NRHS",
                     InfoSpacing, "PadA",
                     InfoSpacing, "PadB",
                     InfoSpacing, "NB",
                     InfoSpacing, "IB",
                     InfoSpacing, "Hous. mode");
        }
        return;
    }
    // Return column values.
    snprintf(info, InfoLen,
             "%*d %*d %*d %*d %*d %*d %*d %*c",
             InfoSpacing, param[PARAM_M].i,
             InfoSpacing, param[PARAM_N].i,
             InfoSpacing, param[PARAM_NRHS].i,
             InfoSpacing, param[PARAM_PADA].i,
             InfoSpacing, param[PARAM_PADB].i,
             InfoSpacing, param[PARAM_NB].i,
             InfoSpacing, param[PARAM_IB].i,
             InfoSpacing, param[PARAM_HMODE].c);

    //================================================================
    // Set parameters.
    //================================================================
    int m    = param[PARAM_M].i;
    int n    = param[PARAM_N].i;
    int nrhs = param[PARAM_NRHS].i;

    int lda = imax(1, m + param[PARAM_PADA].i);
    int ldb = imax(1, imax(m, n) + param[PARAM_PADB].i);

    int test = param[PARAM_TEST].c == 'y';
    float tol = param[PARAM_TOL].d * LAPACKE_slamch('E');

    //================================================================
    // Set tuning parameters.
    //================================================================
    plasma_set(PlasmaNb, param[PARAM_NB].i);
    plasma_set(PlasmaIb, param[PARAM_IB].i);
    if (param[PARAM_HMODE].c == 't') {
        plasma_set(PlasmaHouseholderMode, PlasmaTreeHouseholder);
    }
    else {
        plasma_set(PlasmaHouseholderMode, PlasmaFlatHouseholder);
    }

    //================================================================
    // Allocate and initialize arrays.
    //================================================================
    float *A =
        (float*)malloc((size_t)lda*n*sizeof(float));
    assert(A != NULL);

    float *B =
        (float*)malloc((size_t)ldb*nrhs*
                                    sizeof(float));
    assert(B != NULL);

    int seed[] = {0, 0, 0, 1};
    lapack_int retval;
    retval = LAPACKE_slarnv(1, seed, (size_t)lda*n, A);
    assert(retval == 0);

    retval = LAPACKE_slarnv(1, seed, (size_t)ldb*nrhs, B);
    assert(retval == 0);

    // store the original arrays if residual is to be evaluated
    float *Aref = NULL;
    float *Bref = NULL;
    if (test) {
        Aref = (float*)malloc(
            (size_t)lda*n*sizeof(float));
        assert(Aref != NULL);

        memcpy(Aref, A, (size_t)lda*n*sizeof(float));

        Bref = (float*)malloc(
            (size_t)ldb*nrhs*sizeof(float));
        assert(Bref != NULL);

        memcpy(Bref, B, (size_t)ldb*nrhs*sizeof(float));
    }

    //================================================================
    // Prepare the descriptor for matrix T.
    //================================================================
    plasma_desc_t T;

    //================================================================
    // Run and time PLASMA.
    //================================================================
    plasma_time_t start = omp_get_wtime();
    plasma_sgels(PlasmaNoTrans, m, n, nrhs,
                 A, lda,
                 &T,
                 B, ldb);
    plasma_time_t stop = omp_get_wtime();
    plasma_time_t time = stop-start;

    param[PARAM_TIME].d = time;

    float flop;
    if (m >= n) {
        // cost of QR-based factorization and solve
        flop = flops_sgeqrf(m, n) + flops_sgeqrs(m, n, nrhs);
    }
    else {
        // cost of LQ-based factorization, triangular solve, and Q^T application
        flop = flops_sgelqf(m, n) +
               flops_strsm(PlasmaLeft, m, nrhs) +
               flops_sormlq(PlasmaLeft, n, nrhs, m);
    }
    param[PARAM_GFLOPS].d = flop / time / 1e9;

    //================================================================
    // Test results by solving a linear system.
    //================================================================
    if (test) {
        // |A|_F
        float work[1];
        float Anorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'F', m, n,
                                           Aref, lda, work);

        // |B|_F
        float Bnorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'F', m, nrhs,
                                           Bref, ldb, work);

        // |X|_F, solution X is now stored in the n-by-nrhs part of B
        float Xnorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'F', n, nrhs,
                                           B, ldb, work);

        // compute residual and store it in B = A*X - B
        float zone  =  1.0;
        float zmone = -1.0;
        float zzero =  0.0;
        cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, nrhs, n,
                    (zone), Aref, lda, B, ldb,
                    (zmone), Bref, ldb);

        // Compute B = A^T * (A*X - B)
        cblas_sgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, n, nrhs, m,
                    (zone), Aref, lda, Bref, ldb,
                    (zzero), B, ldb);

        // |RES|_F
        float Rnorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'F', n, nrhs,
                                           B, ldb, work);

        // normalize the result
        float result = Rnorm / ((Anorm*Xnorm+Bnorm)*n);

        param[PARAM_ERROR].d = result;
        param[PARAM_SUCCESS].i = result < tol;
    }

    //================================================================
    // Free arrays.
    //================================================================
    plasma_desc_destroy(&T);
    free(A);
    free(B);
    if (test) {
        free(Aref);
        free(Bref);
    }
}
