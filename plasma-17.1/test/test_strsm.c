/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from test/test_ztrsm.c, normal z -> s, Mon Nov 22 19:22:48 2021
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

#include <omp.h>

#define REAL

#define A(i_, j_) A[(i_) + (size_t)lda*(j_)]

/***************************************************************************//**
 *
 * @brief Tests STRSM.
 *
 * @param[in]  param - array of parameters
 * @param[out] info  - string of column labels or column values; length InfoLen
 *
 * If param is NULL and info is NULL,     print usage and return.
 * If param is NULL and info is non-NULL, set info to column labels and return.
 * If param is non-NULL and info is non-NULL, set info to column values
 * and run test.
 ******************************************************************************/
void test_strsm(param_value_t param[], char *info)
{
    //================================================================
    // Print usage info or return column labels or values.
    //================================================================
    if (param == NULL) {
        if (info == NULL) {
            // Print usage info.
            print_usage(PARAM_SIDE);
            print_usage(PARAM_UPLO);
            print_usage(PARAM_TRANSA);
            print_usage(PARAM_DIAG);
            print_usage(PARAM_M);
            print_usage(PARAM_N);
            print_usage(PARAM_ALPHA);
            print_usage(PARAM_PADA);
            print_usage(PARAM_PADB);
            print_usage(PARAM_NB);
        }
        else {
            // Return column labels.
            snprintf(info, InfoLen,
                     "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s",
                     InfoSpacing, "side",
                     InfoSpacing, "uplo",
                     InfoSpacing, "TransA",
                     InfoSpacing, "diag",
                     InfoSpacing, "M",
                     InfoSpacing, "N",
                     InfoSpacing, "alpha",
                     InfoSpacing, "PadA",
                     InfoSpacing, "PadB",
                     InfoSpacing, "NB");
        }
        return;
    }
    // Return column values.
    snprintf(info, InfoLen,
             "%*c %*c %*c %*c %*d %*d %*.4f %*d %*d %*d",
             InfoSpacing, param[PARAM_SIDE].c,
             InfoSpacing, param[PARAM_UPLO].c,
             InfoSpacing, param[PARAM_TRANSA].c,
             InfoSpacing, param[PARAM_DIAG].c,
             InfoSpacing, param[PARAM_M].i,
             InfoSpacing, param[PARAM_N].i,
             InfoSpacing, creal(param[PARAM_ALPHA].z),
             InfoSpacing, param[PARAM_PADA].i,
             InfoSpacing, param[PARAM_PADB].i,
             InfoSpacing, param[PARAM_NB].i);

    //================================================================
    // Set parameters.
    //================================================================
    plasma_enum_t side = plasma_side_const(param[PARAM_SIDE].c);
    plasma_enum_t uplo = plasma_uplo_const(param[PARAM_UPLO].c);
    plasma_enum_t transa = plasma_trans_const(param[PARAM_TRANSA].c);
    plasma_enum_t diag = plasma_diag_const(param[PARAM_DIAG].c);

    int m = param[PARAM_M].i;
    int n = param[PARAM_N].i;
    int Am;

    if (side == PlasmaLeft) {
        Am = m;
    }
    else {
        Am = n;
    }

    int lda = imax(1, Am + param[PARAM_PADA].i);
    int ldb = imax(1, m + param[PARAM_PADB].i);

    int test = param[PARAM_TEST].c == 'y';
    float tol = param[PARAM_TOL].d * LAPACKE_slamch('E');

#ifdef COMPLEX
    float alpha = param[PARAM_ALPHA].z;
#else
    float alpha = creal(param[PARAM_ALPHA].z);
#endif

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

    float *B =
        (float*)malloc((size_t)ldb*n*sizeof(float));
    assert(B != NULL);

    int *ipiv;
    ipiv = (int*)malloc((size_t)Am*sizeof(int));
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

    LAPACKE_sgetrf(CblasColMajor, Am, Am, A, lda, ipiv);

    if (diag == PlasmaUnit && uplo == PlasmaUpper) {
        // U = L^T
        for (int j = 0; j < Am; j++) {
            for (int i = 0; i < j; i++) {
                A(i,j) = A(j,i);
            }
        }
    }
    else if (diag == PlasmaNonUnit && uplo == PlasmaLower) {
        // L = U^T
        for (int j = 0; j < Am; j++) {
            for (int i = 0; i < j; i++) {
                A(j,i) = A(i,j);
            }
        }
    }

    retval = LAPACKE_slarnv(1, seed, (size_t)ldb*n, B);
    assert(retval == 0);

    float *Bref = NULL;
    if (test) {
        Bref = (float*)malloc(
            (size_t)ldb*n*sizeof(float));
        assert(Bref != NULL);

        memcpy(Bref, B, (size_t)ldb*n*sizeof(float));
    }

    //================================================================
    // Run and time PLASMA.
    //================================================================
    plasma_time_t start = omp_get_wtime();

    plasma_strsm(
        side, uplo,
        transa, diag,
        m, n,
        alpha, A, lda,
               B, ldb);

    plasma_time_t stop = omp_get_wtime();
    plasma_time_t time = stop-start;

    param[PARAM_TIME].d = time;
    param[PARAM_GFLOPS].d = flops_strsm(side, m, n) / time / 1e9;

    //================================================================
    // Test results by checking the residual
    // ||alpha*B - A*X|| / (||A||*||X||)
    //================================================================
    if (test) {
        float zone  =  1.0;
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
                                     &Am, &Am, A, &lda, work);
        //float Anorm = LAPACKE_slantr_work(
        //                   LAPACK_COL_MAJOR, 'F', lapack_const(uplo),
        //                   lapack_const(diag), Am, Am, A, lda, work);

        float Xnorm = LAPACKE_slange_work(
                           LAPACK_COL_MAJOR, 'F', m, n, B, ldb, work);

        // B = A*X
        cblas_strmm(
            CblasColMajor,
            (CBLAS_SIDE)side, (CBLAS_UPLO)uplo,
            (CBLAS_TRANSPOSE)transa, (CBLAS_DIAG)diag,
            m, n,
            (zone), A, lda,
            B, ldb);

        // B -= alpha*Bref
        cblas_sscal((size_t)ldb*n, (alpha), Bref, 1);
        cblas_saxpy((size_t)ldb*n, (zmone), Bref, 1, B, 1);

        float error = LAPACKE_slange_work(
                           LAPACK_COL_MAJOR, 'F', m, n, B, ldb, work);
        if (Anorm * Xnorm != 0)
            error /= (Anorm * Xnorm);

        param[PARAM_ERROR].d = error;
        param[PARAM_SUCCESS].i = error < tol;
    }

    //================================================================
    // Free arrays.
    //================================================================
    free(A);
    free(B);
    free(ipiv);
    if (test)
        free(Bref);
}
