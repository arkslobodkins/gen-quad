/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from test/test_zgelqs.c, normal z -> c, Mon Nov 22 19:22:50 2021
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

#define COMPLEX

/***************************************************************************//**
 *
 * @brief Tests CGELQS.
 *
 * @param[in]  param - array of parameters
 * @param[out] info  - string of column labels or column values; length InfoLen
 *
 * If param is NULL and info is NULL,     print usage and return.
 * If param is NULL and info is non-NULL, set info to column labels and return.
 * If param is non-NULL and info is non-NULL, set info to column values
 * and run test.
 ******************************************************************************/
void test_cgelqs(param_value_t param[], char *info)
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
    plasma_complex32_t *A =
        (plasma_complex32_t*)malloc((size_t)lda*n*sizeof(plasma_complex32_t));
    assert(A != NULL);

    plasma_complex32_t *B =
        (plasma_complex32_t*)malloc((size_t)ldb*nrhs*
                                    sizeof(plasma_complex32_t));
    assert(B != NULL);

    int seed[] = {0, 0, 0, 1};
    lapack_int retval;
    retval = LAPACKE_clarnv(1, seed, (size_t)lda*n, A);
    assert(retval == 0);

    retval = LAPACKE_clarnv(1, seed, (size_t)ldb*nrhs, B);
    assert(retval == 0);

    // store the original arrays if residual is to be evaluated
    plasma_complex32_t *Aref = NULL;
    plasma_complex32_t *Bref = NULL;
    if (test) {
        Aref = (plasma_complex32_t*)malloc(
            (size_t)lda*n*sizeof(plasma_complex32_t));
        assert(Aref != NULL);

        memcpy(Aref, A, (size_t)lda*n*sizeof(plasma_complex32_t));

        Bref = (plasma_complex32_t*)malloc(
            (size_t)ldb*nrhs*sizeof(plasma_complex32_t));
        assert(Bref != NULL);

        memcpy(Bref, B, (size_t)ldb*nrhs*sizeof(plasma_complex32_t));
    }

    //================================================================
    // Prepare the descriptor for matrix T.
    //================================================================
    plasma_desc_t T;

    //================================================================
    // Run and time PLASMA.
    //================================================================
    // prepare LQ factorization of A - only auxiliary for this test,
    // time is not measured
    plasma_cgelqf(m, n, A, lda, &T);

    // perform solution of the system by the prepared LQ factorization of A
    plasma_time_t start = omp_get_wtime();
    plasma_cgelqs(m, n, nrhs,
                  A, lda,
                  T,
                  B, ldb);

    plasma_time_t stop = omp_get_wtime();
    plasma_time_t time = stop-start;

    param[PARAM_TIME].d = time;
    param[PARAM_GFLOPS].d =
        (flops_ctrsm(PlasmaLeft, m, nrhs) +
         flops_cunmlq(PlasmaLeft, n, nrhs, m)) / time / 1e9;

    //================================================================
    // Test results by solving a linear system.
    //================================================================
    if (test) {
        // |A|_F
        float work[1];
        float Anorm = LAPACKE_clange_work(LAPACK_COL_MAJOR, 'F', m, n,
                                           Aref, lda, work);

        // |B|_F
        float Bnorm = LAPACKE_clange_work(LAPACK_COL_MAJOR, 'F', m, nrhs,
                                           Bref, ldb, work);

        // |X|_F, solution X is now stored in the n-by-nrhs part of B
        float Xnorm = LAPACKE_clange_work(LAPACK_COL_MAJOR, 'F', n, nrhs,
                                           B, ldb, work);

        // compute residual and store it in B = A*X - B
        plasma_complex32_t zone  =  1.0;
        plasma_complex32_t zmone = -1.0;
        plasma_complex32_t zzero =  0.0;
        cblas_cgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, nrhs, n,
                    CBLAS_SADDR(zone), Aref, lda, B, ldb,
                    CBLAS_SADDR(zmone), Bref, ldb);

        // Compute B = A^H * (A*X - B)
        cblas_cgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, n, nrhs, m,
                    CBLAS_SADDR(zone), Aref, lda, Bref, ldb,
                    CBLAS_SADDR(zzero), B, ldb);

        // |RES|_F
        float Rnorm = LAPACKE_clange_work(LAPACK_COL_MAJOR, 'F', n, nrhs,
                                           B, ldb, work);

        // normalize the result
        float result = Rnorm / ( (Anorm*Xnorm+Bnorm)*n);

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
