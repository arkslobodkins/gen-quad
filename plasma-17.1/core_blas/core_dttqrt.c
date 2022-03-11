/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from core_blas/core_zttqrt.c, normal z -> d, Thu Mar 10 18:58:54 2022
 *
 **/

#include "core_blas.h"
#include "plasma_types.h"
#include "plasma_internal.h"
#include "core_lapack.h"

#include <omp.h>

// This will be swapped during the automatic code generation.
#undef REAL
#define REAL

/***************************************************************************//**
 *
 * @ingroup core_ttqrt
 *
 * Computes a QR factorization of a rectangular matrix
 * formed by coupling an n-by-n upper triangular tile A1
 * on top of an m-by-n upper triangular tile A2:
 *
 *    | A1 | = Q * R
 *    | A2 |
 *
 *******************************************************************************
 *
 * @param[in] m
 *         The number of columns of the tile A2. m >= 0.
 *
 * @param[in] n
 *         The number of rows of the tile A1.
 *         The number of columns of the tiles A1 and A2. n >= 0.
 *
 * @param[in] ib
 *         The inner-blocking size.  ib >= 0.
 *
 * @param[in,out] A1
 *         On entry, the n-by-n tile A1.
 *         On exit, the elements on and above the diagonal of the array
 *         contain the n-by-n upper trapezoidal tile R;
 *         the elements below the diagonal are not referenced.
 *
 * @param[in] lda1
 *         The leading dimension of the array A1. lda1 >= max(1,n).
 *
 * @param[in,out] A2
 *         On entry, the m-by-n upper triangular tile A2.
 *         On exit, the elements on and above the diagonal of the array
 *         with the matrix T represent
 *         the orthogonal tile Q as a product of elementary reflectors
 *
 * @param[in] lda2
 *         The leading dimension of the tile A2. lda2 >= max(1,m).
 *
 * @param[out] T
 *         The ib-by-n triangular factor T of the block reflector.
 *         T is upper triangular by block (economic storage);
 *         The rest of the array is not referenced.
 *
 * @param[in] ldt
 *         The leading dimension of the array T. ldt >= ib.
 *
 * @param tau
 *         Auxiliary workspace array of length n.
 *
 * @param work
 *         Auxiliary workspace array of length ib*n.
 *
 *******************************************************************************
 *
 * @retval PlasmaSuccess successful exit
 * @retval < 0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
int core_dttqrt(int m, int n, int ib,
                double *A1, int lda1,
                double *A2, int lda2,
                double *T,  int ldt,
                double *tau,
                double *work)
{
    // Check input arguments.
    if (m < 0) {
        coreblas_error("illegal value of m");
        return -1;
    }
    if (n < 0) {
        coreblas_error("illegal value of n");
        return -2;
    }
    if (ib < 0) {
        coreblas_error("illegal value of ib");
        return -3;
    }
    if (A1 == NULL) {
        coreblas_error("NULL A1");
        return -4;
    }
    if (lda1 < imax(1, m) && m > 0) {
        coreblas_error("illegal value of lda1");
        return -5;
    }
    if (A2 == NULL) {
        coreblas_error("NULL A2");
        return -6;
    }
    if (lda2 < imax(1, m) && m > 0) {
        coreblas_error("illegal value of lda2");
        return -7;
    }
    if (T == NULL) {
        coreblas_error("NULL T");
        return -8;
    }
    if (ldt < imax(1, ib) && ib > 0) {
        coreblas_error("illegal value of ldt");
        return -9;
    }
    if (tau == NULL) {
        coreblas_error("NULL tau");
        return -10;
    }
    if (work == NULL) {
        coreblas_error("NULL work");
        return -11;
    }

    // quick return
    if ((m == 0) || (n == 0) || (ib == 0))
        return PlasmaSuccess;

    // TODO: Need to check why some cases require this to avoid
    // uninitialized values
    //core_dlaset(PlasmaGeneral, ib, n, 0.0, 0.0, T, ldt);

    for (int ii = 0; ii < n; ii += ib) {
        int sb = imin(n-ii, ib);
        for (int i = 0; i < sb; i++) {
            int j  = ii + i;
            int mi = imin(j+1, m);
            int ni = sb-i-1;

            // Generate elementary reflector H(j) to annihilate A2(1:mi, j).
            LAPACKE_dlarfg_work(
                    mi+1, &A1[lda1*j+j], &A2[lda2*j], 1, &tau[j]);

            if (ni > 0) {
                // Apply H(j-1) to A(j:m, j+1:ii+ib) from the left.
                cblas_dcopy(
                    ni,
                    &A1[lda1*(j+1)+j], lda1,
                    work, 1);

#ifdef COMPLEX
                LAPACKE_dlacgv_work(ni, work, 1);
#endif
                double zone  = 1.0;
                cblas_dgemv(
                    CblasColMajor, (CBLAS_TRANSPOSE)PlasmaTrans,
                    mi, ni,
                    (zone), &A2[lda2*(j+1)], lda2,
                                       &A2[lda2*j],     1,
                    (zone), work,            1);
#ifdef COMPLEX
                LAPACKE_dlacgv_work(ni, work, 1);
#endif
                double alpha = -(tau[j]);
                cblas_daxpy(
                    ni, (alpha),
                    work, 1,
                    &A1[lda1*(j+1)+j], lda1);
#ifdef COMPLEX
                LAPACKE_dlacgv_work(ni, work, 1);
#endif
                cblas_dger(
                    CblasColMajor, mi, ni,
                    (alpha), &A2[lda2*j], 1,
                    work, 1,
                    &A2[lda2*(j+1)], lda2);
            }

            // Calculate T.
            // T(0:i-1, j) = alpha * A2(0:m-1, ii:j-1)^T * A2(0:m-1, j)
            if (i > 0) {
                int l = imin(i, imax(0, m-ii));
                double alpha = -(tau[j]);

                core_dpemv(
                        PlasmaTrans, PlasmaColumnwise,
                        imin(j, m), i, l,
                        alpha, &A2[lda2*ii], lda2,
                               &A2[lda2*j],  1,
                        0.0, &T[ldt*j],    1,
                        work);

                // T(0:i-1, j) = T(0:i-1, ii:j-1) * T(0:i-1, j)
                cblas_dtrmv(
                        CblasColMajor, (CBLAS_UPLO)PlasmaUpper,
                        (CBLAS_TRANSPOSE)PlasmaNoTrans,
                        (CBLAS_DIAG)PlasmaNonUnit,
                        i, &T[ldt*ii], ldt,
                           &T[ldt*j], 1);
            }
            T[ldt*j+i] = tau[j];
        }

        // Apply Q^T to the rest of the matrix from the left.
        if (n > ii+sb) {
            int mi = imin(ii+sb, m);
            int ni = n-(ii+sb);
            int l  = imin(sb, imax(0, mi-ii));
            core_dparfb(
                PlasmaLeft, PlasmaTrans,
                PlasmaForward, PlasmaColumnwise,
                ib, ni, mi, ni, sb, l,             //replaced sb by ib
                &A1[lda1*(ii+sb)+ii], lda1,
                &A2[lda2*(ii+sb)], lda2,
                &A2[lda2*ii], lda2,
                &T[ldt*ii], ldt,
                work, sb);
        }
    }

    return PlasmaSuccess;
}

/******************************************************************************/
void core_omp_dttqrt(int m, int n, int ib,
                     double *A1, int lda1,
                     double *A2, int lda2,
                     double *T,  int ldt,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request)
{
    #pragma omp task depend(inout:A1[0:lda1*n]) \
                     depend(inout:A2[0:lda2*n]) \
                     depend(out:T[0:ib*n])
    {
        if (sequence->status == PlasmaSuccess) {
            // Prepare workspaces.
            int tid = omp_get_thread_num();
            double *tau = ((double*)work.spaces[tid]);

            // Call the kernel.
            int info = core_dttqrt(m, n, ib,
                                   A1, lda1,
                                   A2, lda2,
                                   T,  ldt,
                                   tau,
                                   tau+n);

            if (info != PlasmaSuccess) {
                plasma_error("core_dttqrt() failed");
                plasma_request_fail(sequence, request, PlasmaErrorInternal);
            }
        }
    }
}
