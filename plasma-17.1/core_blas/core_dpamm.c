/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from core_blas/core_zpamm.c, normal z -> d, Tue Feb  8 19:16:28 2022
 *
 **/

#include "core_blas.h"
#include "plasma_types.h"
#include "plasma_internal.h"
#include "core_lapack.h"

static inline int core_dpamm_a2(plasma_enum_t side, plasma_enum_t trans,
                                plasma_enum_t uplo,
                                int m, int n, int k, int l, int vi2, int vi3,
                                      double *A2, int lda2,
                                const double *V,  int ldv,
                                      double *W,  int ldw);

static inline int core_dpamm_w(plasma_enum_t side, plasma_enum_t trans,
                               plasma_enum_t uplo,
                               int m, int n, int k, int l, int vi2, int vi3,
                               const double *A1, int lda1,
                                     double *A2, int lda2,
                               const double *V,  int ldv,
                                     double *W,  int ldw);

/***************************************************************************//**
 *
 * @ingroup core_pamm
 *
 *  Performs one of the matrix-matrix operations
 *
 *                    PlasmaLeft                PlasmaRight
 *     OP PlasmaW  :  W  = A1 + op(V) * A2  or  W  = A1 + A2 * op(V)
 *     OP PlasmaA2 :  A2 = A2 - op(V) * W   or  A2 = A2 - W * op(V)
 *
 *  where  op( V ) is one of
 *
 *     op( V ) = V   or   op( V ) = V^T
 *
 *  A1, A2 and W are general matrices, and V is:
 *
 *        l = k: rectangle + triangle
 *        l < k: rectangle + trapezoid
 *        l = 0: rectangle
 *
 *  Size of V, both rowwise and columnwise, is:
 *
 *         ----------------------
 *          side   trans    size
 *         ----------------------
 *          left     N     M x K
 *                   T     K x M
 *          right    N     K x N
 *                   T     N x K
 *         ----------------------
 *
 *  PlasmaLeft (columnwise and rowwise):
 *
 *              |    K    |                 |         M         |
 *           _  __________   _              _______________        _
 *              |    |    |                 |             | \
 *     V:       |    |    |           V^T:  |_____________|___\    K
 *              |    |    | M-L             |                  |
 *           M  |    |    |                 |__________________|   _
 *              |____|    |  _
 *              \    |    |                 |    M - L    | L  |
 *                \  |    |  L
 *           _      \|____|  _
 *
 *
 *  PlasmaRight (columnwise and rowwise):
 *
 *          |         K         |                   |    N    |
 *          _______________        _             _  __________   _
 *          |             | \                       |    |    |
 *    V^T:  |_____________|___\    N        V:      |    |    |
 *          |                  |                    |    |    | K-L
 *          |__________________|   _             K  |    |    |
 *                                                  |____|    |  _
 *          |    K - L    | L  |                    \    |    |
 *                                                    \  |    |  L
 *                                               _      \|____|  _
 *
 *  Arguments
 *  ==========
 *
 * @param[in] op
 *
 *         OP specifies which operation to perform:
 *
 *         - PlasmaW  : W  = A1 + op(V) * A2  or  W  = A1 + A2 * op(V)
 *         - PlasmaA2 : A2 = A2 - op(V) * W   or  A2 = A2 - W * op(V)
 *
 * @param[in] side
 *
 *         SIDE specifies whether  op( V ) multiplies A2
 *         or W from the left or right as follows:
 *
 *         - PlasmaLeft  : multiply op( V ) from the left
 *                            OP PlasmaW  :  W  = A1 + op(V) * A2
 *                            OP PlasmaA2 :  A2 = A2 - op(V) * W
 *
 *         - PlasmaRight : multiply op( V ) from the right
 *                            OP PlasmaW  :  W  = A1 + A2 * op(V)
 *                            OP PlasmaA2 :  A2 = A2 - W * op(V)
 *
 * @param[in] storev
 *
 *         Indicates how the vectors which define the elementary
 *         reflectors are stored in V:
 *
 *         - PlasmaColumnwise
 *         - PlasmaRowwise
 *
 * @param[in] m
 *         The number of rows of the A1, A2 and W
 *         If SIDE is PlasmaLeft, the number of rows of op( V )
 *
 * @param[in] n
 *         The number of columns of the A1, A2 and W
 *         If SIDE is PlasmaRight, the number of columns of op( V )
 *
 * @param[in] k
 *         If SIDE is PlasmaLeft, the number of columns of op( V )
 *         If SIDE is PlasmaRight, the number of rows of op( V )
 *
 * @param[in] l
 *         The size of the triangular part of V
 *
 * @param[in] A1
 *         On entry, the m-by-n tile A1.
 *
 * @param[in] lda1
 *         The leading dimension of the array A1. lda1 >= max(1,m).
 *
 * @param[in,out] A2
 *         On entry, the m-by-n tile A2.
 *         On exit, if OP is PlasmaA2 A2 is overwritten
 *
 * @param[in] lda2
 *         The leading dimension of the tile A2. lda2 >= max(1,m).
 *
 * @param[in] V
 *         The matrix V as described above.
 *         If SIDE is PlasmaLeft : op( V ) is m-by-k
 *         If SIDE is PlasmaRight: op( V ) is k-by-n
 *
 * @param[in] ldv
 *         The leading dimension of the array V.
 *
 * @param[in,out] W
 *         On entry, the m-by-n matrix W.
 *         On exit, W is overwritten either if OP is PlasmaA2 or PlasmaW.
 *         If OP is PlasmaA2, W is an input and is used as a workspace.
 *
 * @param[in] ldw
 *         The leading dimension of array W.
 *
 *******************************************************************************
 *
 * @retval PlasmaSuccess successful exit
 * @retval < 0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
int core_dpamm(plasma_enum_t op, plasma_enum_t side, plasma_enum_t storev,
               int m, int n, int k, int l,
               const double *A1, int lda1,
                     double *A2, int lda2,
               const double *V,  int ldv,
                     double *W,  int ldw)
{
    // Check input arguments.
    if ((op != PlasmaW) && (op != PlasmaA2)) {
        coreblas_error("illegal value of op");
        return -1;
    }
    if ((side != PlasmaLeft) && (side != PlasmaRight)) {
        coreblas_error("illegal value of side");
        return -2;
    }
    if ((storev != PlasmaColumnwise) && (storev != PlasmaRowwise)) {
        coreblas_error("illegal value of storev");
        return -3;
    }
    if (m < 0) {
        coreblas_error("illegal value of m");
        return -4;
    }
    if (n < 0) {
        coreblas_error("illegal value of n");
        return -5;
    }
    if (k < 0) {
        coreblas_error("illegal value of k");
        return -6;
    }
    if (l < 0) {
        coreblas_error("illegal value of l");
        return -7;
    }
    if (A1 == NULL) {
        coreblas_error("NULL A1");
        return -8;
    }
    if (lda1 < 0) {
        coreblas_error("illegal value of lda1");
        return -9;
    }
    if (A2 == NULL) {
        coreblas_error("NULL A2");
        return -10;
    }
    if (lda2 < 0) {
        coreblas_error("illegal value of lda2");
        return -11;
    }
    if (V == NULL) {
        coreblas_error("NULL V");
        return -12;
    }
    if (ldv < 0) {
        coreblas_error("illegal value of ldv");
        return -13;
    }
    if (W == NULL) {
        coreblas_error("NULL W");
        return -14;
    }
    if (ldw < 0) {
        coreblas_error("illegal value of ldw");
        return -15;
    }

    // quick return
    if (m == 0 || n == 0 || k == 0)
        return PlasmaSuccess;

    // TRANS is set as:
    //
    //        -------------------------------------
    //         side   direct     PlasmaW  PlasmaA2
    //        -------------------------------------
    //         left   colwise       T        N
    //                rowwise       N        T
    //         right  colwise       N        T
    //                rowwise       T        N
    //        -------------------------------------

    plasma_enum_t uplo;
    int trans;
    int vi2, vi3;

    //===================
    // PlasmaColumnwise
    //===================
    if (storev == PlasmaColumnwise) {
        uplo = CblasUpper;
        if (side == PlasmaLeft) {
            trans = op == PlasmaA2 ? PlasmaNoTrans : PlasmaTrans;
            vi2 = trans == PlasmaNoTrans ? m - l : k - l;
        }
        else {
            trans = op == PlasmaW ? PlasmaNoTrans : PlasmaTrans;
            vi2 = trans == PlasmaNoTrans ? k - l : n - l;
        }
        vi3 = ldv * l;
    }
    //================
    // PlasmaRowwise
    //================
    else {
        uplo = CblasLower;
        if (side == PlasmaLeft) {
            trans = op == PlasmaW ? PlasmaNoTrans : PlasmaTrans;
            vi2 = trans == PlasmaNoTrans ? k - l : m - l;
        }
        else {
            trans = op == PlasmaA2 ? PlasmaNoTrans : PlasmaTrans;
            vi2 = trans == PlasmaNoTrans ? n - l : k - l;
        }
        vi2 *= ldv;
        vi3  = l;
    }

    if (op == PlasmaW) {
        core_dpamm_w(side, trans, uplo,
                     m, n, k, l, vi2, vi3,
                     A1, lda1,
                     A2, lda2,
                     V, ldv,
                     W, ldw);
    }
    else if (op == PlasmaA2) {
        core_dpamm_a2(side, trans, uplo,
                      m, n, k, l, vi2, vi3,
                      A2, lda2,
                      V,  ldv,
                      W,  ldw);
    }

    return PlasmaSuccess;
}

/******************************************************************************/
static inline int core_dpamm_w(
        plasma_enum_t side, plasma_enum_t trans, plasma_enum_t uplo,
        int m, int n, int k, int l, int vi2, int vi3,
        const double *A1, int lda1,
              double *A2, int lda2,
        const double *V,  int ldv,
              double *W,  int ldw)
{
    // W = A1 + op(V) * A2  or  W = A1 + A2 * op(V)

    double zone  = 1.0;
    double zzero = 0.0;

    //=============
    // PlasmaLeft
    //=============
    if (side == PlasmaLeft) {
        if ((trans == PlasmaTrans && uplo == CblasUpper) ||
            (trans == PlasmaNoTrans && uplo == CblasLower)) {
            // W = A1 + V^T * A2

            // W = A2_2
            LAPACKE_dlacpy_work(LAPACK_COL_MAJOR,
                                lapack_const(PlasmaGeneral),
                                l, n,
                                &A2[k-l], lda2,
                                 W,       ldw);

            // W = V_2^T * W + V_1^T * A2_1 (ge+tr, top L rows of V^T)
            if (l > 0) {
                // W = V_2^T * W
                cblas_dtrmm(CblasColMajor,
                            CblasLeft, (CBLAS_UPLO)uplo,
                            (CBLAS_TRANSPOSE)trans, CblasNonUnit,
                            l, n,
                            (zone), &V[vi2], ldv,
                                                W,      ldw);

                // W = W + V_1^T * A2_1
                if (k > l) {
                    cblas_dgemm(CblasColMajor,
                                (CBLAS_TRANSPOSE)trans, CblasNoTrans,
                                l, n, k-l,
                                (zone), V,  ldv,
                                                   A2, lda2,
                                (zone), W,  ldw);
                }
            }

            // W_2 = V_3^T * A2: (ge, bottom M-L rows of V^T)
            if (m > l) {
                cblas_dgemm(CblasColMajor,
                            (CBLAS_TRANSPOSE)trans, CblasNoTrans,
                            (m-l), n, k,
                            (zone),  &V[vi3], ldv,
                                                 A2,     lda2,
                            (zzero), &W[l],   ldw);
            }

            // W = A1 + W
            for (int j = 0; j < n; j++) {
                cblas_daxpy(m, (zone),
                            &A1[lda1*j], 1,
                            &W[ldw*j], 1);
            }
        }
        else {
            coreblas_error(
                "Left Upper/NoTrans & Lower/[Conj]Trans not implemented");
            return PlasmaErrorNotSupported;
        }
    }
    //==============
    // PlasmaRight
    //==============
    else {
        if ((trans == PlasmaTrans && uplo == CblasUpper) ||
            (trans == PlasmaNoTrans && uplo == CblasLower)) {
            coreblas_error(
                "Right Upper/[Conj]Trans & Lower/NoTrans not implemented");
            return PlasmaErrorNotSupported;
        }
        else {
            // W = A1 + A2 * V
            if (l > 0) {
                // W = A2_2
                LAPACKE_dlacpy_work(LAPACK_COL_MAJOR,
                                    lapack_const(PlasmaGeneral),
                                    m, l,
                                    &A2[lda2*(k-l)], lda2,
                                     W,              ldw);

                // W = W * V_2 --> W = A2_2 * V_2
                cblas_dtrmm(CblasColMajor,
                            CblasRight, (CBLAS_UPLO)uplo,
                            (CBLAS_TRANSPOSE)trans, CblasNonUnit,
                            m, l,
                            (zone), &V[vi2], ldv,
                                                W,      ldw);

                // W = W + A2_1 * V_1
                if (k > l) {
                    cblas_dgemm(CblasColMajor,
                                CblasNoTrans, (CBLAS_TRANSPOSE)trans,
                                m, l, k-l,
                                (zone), A2, lda2,
                                                   V,  ldv,
                                (zone), W,  ldw);
                }
            }

            // W = W + A2 * V_3
            if (n > l) {
                cblas_dgemm(CblasColMajor,
                            CblasNoTrans, (CBLAS_TRANSPOSE)trans,
                            m, n-l, k,
                            (zone),   A2,       lda2,
                                                &V[vi3],   ldv,
                            (zzero), &W[ldw*l], ldw);
            }

            // W = A1 + W
            for (int j = 0; j < n; j++) {
                cblas_daxpy(m, (zone),
                            &A1[lda1*j], 1,
                            &W[ldw*j],   1);
            }
        }
    }

    return PlasmaSuccess;
}

/******************************************************************************/
static inline int core_dpamm_a2(
        plasma_enum_t side, plasma_enum_t trans, plasma_enum_t uplo,
        int m, int n, int k, int l, int vi2, int vi3,
              double *A2, int lda2,
        const double *V,  int ldv,
              double *W,  int ldw)
{
    // A2 = A2 + op(V) * W  or  A2 = A2 + W * op(V)

    double zone  =  1.0;
    double zmone = -1.0;

    //=============
    // PlasmaLeft
    //=============
    if (side == PlasmaLeft) {
        if ((trans == PlasmaTrans && uplo == CblasUpper) ||
            (trans == PlasmaNoTrans && uplo == CblasLower)) {
            coreblas_error(
                "Left Upper/[Conj]Trans & Lower/NoTrans not implemented");
            return PlasmaErrorNotSupported;
        }
        else {
            // A2 = A2 - V * W

            // A2_1 = A2_1 - V_1  * W_1
            if (m > l) {
                cblas_dgemm(CblasColMajor,
                            (CBLAS_TRANSPOSE)trans, CblasNoTrans,
                            m-l, n, l,
                            (zmone), V,  ldv,
                                                W,  ldw,
                            (zone),  A2, lda2);
            }

            // W_1 = V_2 * W_1
            cblas_dtrmm(CblasColMajor,
                        CblasLeft, (CBLAS_UPLO)uplo,
                        (CBLAS_TRANSPOSE)trans, CblasNonUnit,
                        l, n,
                        (zone), &V[vi2], ldv,
                                            W,      ldw);

            // A2_2 = A2_2 - W_1
            for (int j = 0; j < n; j++) {
                cblas_daxpy(l, (zmone),
                            &W[ldw*j], 1,
                            &A2[lda2*j+(m-l)], 1);
            }

            // A2 = A2 - V_3  * W_2
            if (k > l) {
                cblas_dgemm(CblasColMajor,
                            (CBLAS_TRANSPOSE)trans, CblasNoTrans,
                            m, n, (k-l),
                            (zmone), &V[vi3], ldv,
                                                &W[l],   ldw,
                            (zone),   A2,     lda2);
            }
        }
    }
    //==============
    // PlasmaRight
    //==============
    else {
        if ((trans == PlasmaTrans && uplo == CblasUpper) ||
            (trans == PlasmaNoTrans && uplo == CblasLower)) {
            // A2 = A2 - W * V^T

            // A2 = A2 - W_2 * V_3^T
            if (k > l) {
                cblas_dgemm(CblasColMajor,
                            CblasNoTrans, (CBLAS_TRANSPOSE)trans,
                            m, n, k-l,
                            (zmone), &W[ldw*l], ldw,
                                                &V[vi3],   ldv,
                            (zone),   A2,       lda2);
            }

            // A2_1 = A2_1 - W_1 * V_1^T
            if (n > l) {
                cblas_dgemm(CblasColMajor,
                            CblasNoTrans, (CBLAS_TRANSPOSE)trans,
                            m, n-l, l,
                            (zmone), W,  ldw,
                                                V,  ldv,
                            (zone),  A2, lda2);
            }

            // A2_2 =  A2_2 -  W_1 * V_2^T
            if (l > 0) {
                cblas_dtrmm(CblasColMajor,
                            CblasRight, (CBLAS_UPLO)uplo,
                            (CBLAS_TRANSPOSE)trans, CblasNonUnit,
                            m, l,
                            (zmone), &V[vi2], ldv,
                                                 W,      ldw);

                for (int j = 0; j < l; j++) {
                    cblas_daxpy(m, (zone),
                                &W[ldw*j], 1,
                                &A2[lda2*(n-l+j)], 1);
                }
            }
        }
        else {
            coreblas_error(
                "Right Upper/NoTrans & Lower/[Conj]Trans not implemented");
            return PlasmaErrorNotSupported;
        }
    }

    return PlasmaSuccess;
}
