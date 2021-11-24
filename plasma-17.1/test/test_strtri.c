/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from test/test_ztrtri.c, normal z -> s, Mon Nov 22 19:22:54 2021
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
void test_strtri(param_value_t param[], char *info)
{
    //================================================================
    // Print usage info or return column labels or values.
    //================================================================
    if (param == NULL) {
        if (info == NULL) {
            // Print usage info.
            print_usage(PARAM_UPLO);
            print_usage(PARAM_DIAG);
            print_usage(PARAM_N);
            print_usage(PARAM_PADA);
            print_usage(PARAM_NB);
        }
        else {
            // Return column labels.
            snprintf(info, InfoLen,
                     "%*s %*s %*s %*s %*s",
                     InfoSpacing, "uplo",
                     InfoSpacing, "diag",
                     InfoSpacing, "N",
                     InfoSpacing, "PadA",
                     InfoSpacing, "NB");
        }
        return;
    }
    // Return column values.
    snprintf(info, InfoLen,
             "%*c %*c %*d %*d %*d",
             InfoSpacing, param[PARAM_UPLO].c,
             InfoSpacing, param[PARAM_DIAG].c,
             InfoSpacing, param[PARAM_N].i,
             InfoSpacing, param[PARAM_PADA].i,
             InfoSpacing, param[PARAM_NB].i);

    //================================================================
    // Set parameters.
    //================================================================
    plasma_enum_t uplo = plasma_uplo_const(param[PARAM_UPLO].c);
    plasma_enum_t diag = plasma_diag_const(param[PARAM_DIAG].c);

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

    float *Aref;

    int *ipiv;
    ipiv = (int*)malloc((size_t)lda*sizeof(int));
    assert(ipiv != NULL);

    int seed[] = {0, 0, 0, 1};
    lapack_int retval;

    //=================================================================
    // Initialize the matrices.
    // Factor A into LU to get well-conditioned triangular matrices.
    // Use L for unit triangle, and U for non-unit triangle,
    // transposing as necessary.
    // (There is some danger, as L^T or U^T may be much worse conditioned
    // than L or U, but in practice it seems okay.
    // See Higham, Accuracy and Stability of Numerical Algorithms, ch 8.)
    //=================================================================
    retval = LAPACKE_slarnv(1, seed, (size_t)lda*n, A);
    assert(retval == 0);

    LAPACKE_sgetrf(CblasColMajor, n, n, A, lda, ipiv);

    if (diag == PlasmaUnit && uplo == PlasmaUpper) {
        // U = L^T
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < j; i++) {
                A(i, j) = A(j, i);
            }
        }
    }
    else if (diag == PlasmaNonUnit && uplo == PlasmaLower) {
        // L = U^T
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < j; i++) {
                A(j, i) = A(i, j);
            }
        }
    }

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

    plasma_strtri(
        uplo, diag,
        n, A, lda);

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
        float work[1];

        // LAPACKE_[ds]lantr_work has a bug (returns 0)
        // in MKL <= 11.3.3 (at least). Fixed in LAPACK 3.6.1.
        // For now, call LAPACK directly.
        // LAPACK_slantr is a macro for correct name mangling (e.g.
        // adding _ at the end) of the Fortran symbol.
        // The macro is either defined in lapacke.h, or in the file
        // core_lapack_s.h for the use with MKL.
        char normc = 'F';
        char uploc = lapack_const(uplo);
        char diagc = lapack_const(diag);
        float Anorm = LAPACK_slantr(&normc, &uploc, &diagc,
                                     &n, &n, Aref, &lda, work);

        // B = inv(A)
        LAPACKE_strtri_work(
            CblasColMajor,
            lapack_const(uplo), lapack_const(diag),
            n, Aref, lda);

        float Inorm = LAPACK_slantr(&normc, &uploc, &diagc,
                                     &n, &n, Aref, &lda, work);

        // A <- A - Aref
        cblas_saxpy((size_t)lda*n, (zmone), Aref, 1, A, 1);

        char ndiagc = 'N'; // A-Aref cannot have a unit diagonal
        float error = LAPACK_slantr(&normc, &uploc, &ndiagc,
                                     &n, &n, A, &lda, work);
        if (Anorm*Inorm != 0.0)
            error /= (Anorm * Inorm);

        param[PARAM_ERROR].d = error;
        param[PARAM_SUCCESS].i = error < tol;
    }

    //================================================================
    // Free arrays.
    //================================================================
    free(A);
    free(ipiv);
    if (test)
        free(Aref);
}
