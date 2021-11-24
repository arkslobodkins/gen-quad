/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from test/test_zgelqf.c, normal z -> c, Mon Nov 22 19:23:07 2021
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
 * @brief Tests CGELQF.
 *
 * @param[in]  param - array of parameters
 * @param[out] info  - string of column labels or column values; length InfoLen
 *
 * If param is NULL and info is NULL,     print usage and return.
 * If param is NULL and info is non-NULL, set info to column labels and return.
 * If param is non-NULL and info is non-NULL, set info to column values
 * and run test.
 ******************************************************************************/
void test_cgelqf(param_value_t param[], char *info)
{
    //================================================================
    // Print usage info or return column labels or values.
    //================================================================
    if (param == NULL) {
        if (info == NULL) {
            // Print usage info.
            print_usage(PARAM_M);
            print_usage(PARAM_N);
            print_usage(PARAM_PADA);
            print_usage(PARAM_NB);
            print_usage(PARAM_IB);
            print_usage(PARAM_HMODE);
        }
        else {
            // Return column labels.
            snprintf(info, InfoLen,
                     "%*s %*s %*s %*s %*s %*s %*s",
                     InfoSpacing, "M",
                     InfoSpacing, "N",
                     InfoSpacing, "PadA",
                     InfoSpacing, "NB",
                     InfoSpacing, "IB",
                     InfoSpacing, "Hous. mode",
                     InfoSpacing, "Ortho.");
        }
        return;
    }
    // Return column values.
    // ortho. column appended later.
    snprintf(info, InfoLen,
             "%*d %*d %*d %*d %*d %*c",
             InfoSpacing, param[PARAM_M].i,
             InfoSpacing, param[PARAM_N].i,
             InfoSpacing, param[PARAM_PADA].i,
             InfoSpacing, param[PARAM_NB].i,
             InfoSpacing, param[PARAM_IB].i,
             InfoSpacing, param[PARAM_HMODE].c);

    //================================================================
    // Set parameters.
    //================================================================
    int m = param[PARAM_M].i;
    int n = param[PARAM_N].i;

    int lda = imax(1, m + param[PARAM_PADA].i);

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

    int seed[] = {0, 0, 0, 1};
    lapack_int retval;
    retval = LAPACKE_clarnv(1, seed, (size_t)lda*n, A);
    assert(retval == 0);

    plasma_complex32_t *Aref = NULL;
    if (test) {
        Aref = (plasma_complex32_t*)malloc(
            (size_t)lda*n*sizeof(plasma_complex32_t));
        assert(Aref != NULL);

        memcpy(Aref, A, (size_t)lda*n*sizeof(plasma_complex32_t));
    }

    //================================================================
    // Prepare the descriptor for matrix T.
    //================================================================
    plasma_desc_t T;

    //================================================================
    // Run and time PLASMA.
    //================================================================
    plasma_time_t start = omp_get_wtime();
    plasma_cgelqf(m, n, A, lda, &T);
    plasma_time_t stop = omp_get_wtime();
    plasma_time_t time = stop-start;

    param[PARAM_TIME].d = time;
    param[PARAM_GFLOPS].d = flops_cgelqf(m, n) / time / 1e9;

    //=================================================================
    // Test results by checking orthogonality of Q and precision of L*Q
    //=================================================================
    if (test) {
        // Check the orthogonality of Q.
        int minmn = imin(m, n);

        // Allocate space for Q.
        int ldq = minmn;
        plasma_complex32_t *Q =
            (plasma_complex32_t *)malloc((size_t)ldq*n*
                                         sizeof(plasma_complex32_t));
        assert(Q != NULL);

        // Build Q.
        plasma_cunglq(minmn, n, minmn, A, lda, T, Q, ldq);

        // Build the identity matrix
        plasma_complex32_t *Id =
            (plasma_complex32_t *) malloc((size_t)minmn*minmn*
                                          sizeof(plasma_complex32_t));
        assert(Id != NULL);
        LAPACKE_claset_work(LAPACK_COL_MAJOR, 'g', minmn, minmn,
                            0.0, 1.0, Id, minmn);

        // Perform Id - Q * Q^H
        cblas_cherk(CblasColMajor, CblasUpper, CblasNoTrans, minmn, n,
                    -1.0, Q, ldq, 1.0, Id, minmn);

        // work array of size m is needed for computing L_oo norm
        float *work = (float *) malloc((size_t)m*sizeof(float));
        assert(work != NULL);

        // |Id - Q * Q^H|_oo
        float ortho = LAPACKE_clanhe_work(LAPACK_COL_MAJOR, 'I', 'u',
                                           minmn, Id, minmn, work);

        // normalize the result
        // |Id - Q * Q^H|_oo / n
        ortho /= minmn;

        free(Q);
        free(Id);

        // Check the accuracy of A - L * Q
        // LAPACK version does not construct Q, it uses only application of it

        // Extract the L.
        plasma_complex32_t *L =
            (plasma_complex32_t *)malloc((size_t)m*n*
                                         sizeof(plasma_complex32_t));
        assert(L != NULL);
        LAPACKE_claset_work(LAPACK_COL_MAJOR, 'u', m, n,
                            0.0, 0.0, L, m);
        LAPACKE_clacpy_work(LAPACK_COL_MAJOR, 'l', m, n, A, lda, L, m);

        // Compute L * Q.
        plasma_cunmlq(PlasmaRight, PlasmaNoTrans, m, n, minmn, A, lda, T,
                      L, m);

        // Compute the difference.
        // L = A - L*Q
        for (int j = 0; j < n; j++)
            for (int i = 0; i < m; i++)
                L[j*m+i] = Aref[j*lda+i] - L[j*m+i];

        // |A|_oo
        float normA = LAPACKE_clange_work(LAPACK_COL_MAJOR, 'I', m, n,
                                           Aref, lda, work);

        // |A - L*Q|_oo
        float error = LAPACKE_clange_work(LAPACK_COL_MAJOR, 'I', m, n,
                                           L, m, work);

        // normalize the result
        // |A-LQ|_oo / (|A|_oo * n)
        error /= (normA * n);

        param[PARAM_ERROR].d = error;
        param[PARAM_ORTHO].d = ortho;
        param[PARAM_SUCCESS].i = (error < tol && ortho < tol);

        free(work);
        free(L);

        // Return ortho. column value.
        int len = strlen(info);
        snprintf(&info[len], imax(0, InfoLen - len),
                 " %*.2e",
                 InfoSpacing, param[PARAM_ORTHO].d);
    }
    else {
        // No ortho. test.
        int len = strlen(info);
        snprintf(&info[len], imax(0, InfoLen - len),
                 " %*s",
                 InfoSpacing, "---");
    }

    //================================================================
    // Free arrays.
    //================================================================
    plasma_desc_destroy(&T);
    free(A);
    if (test)
        free(Aref);
}
