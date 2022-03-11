/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from compute/pzsymm.c, normal z -> d, Thu Mar 10 18:57:57 2022
 *
 **/

#include "plasma_async.h"
#include "plasma_context.h"
#include "plasma_descriptor.h"
#include "plasma_internal.h"
#include "plasma_types.h"
#include "plasma_workspace.h"
#include "core_blas.h"

#define A(m, n) (double*)plasma_tile_addr(A, m, n)
#define B(m, n) (double*)plasma_tile_addr(B, m, n)
#define C(m, n) (double*)plasma_tile_addr(C, m, n)

/***************************************************************************//**
 *  Parallel tile symmetric matrix-matrix multiplication.
 *  @see plasma_omp_dsymm
 ******************************************************************************/
void plasma_pdsymm(plasma_enum_t side, plasma_enum_t uplo,
                   double alpha, plasma_desc_t A,
                                             plasma_desc_t B,
                   double beta,  plasma_desc_t C,
                   plasma_sequence_t *sequence, plasma_request_t *request)
{
    // Return if failed sequence.
    if (sequence->status != PlasmaSuccess)
        return;

    for (int m = 0; m < C.mt; m++) {
        int mvcm = plasma_tile_mview(C, m);
        int ldcm = plasma_tile_mmain(C, m);
        for (int n = 0; n < C.nt; n++) {
            int nvcn = plasma_tile_nview(C, n);
            if (side == PlasmaLeft) {
                int ldam = plasma_tile_mmain(A, m);
                //===========================
                // PlasmaLeft / PlasmaLower
                //===========================
                if (uplo == PlasmaLower) {
                    for (int k = 0; k < C.mt; k++) {
                        int mvck = plasma_tile_mview(C, k);
                        int ldak = plasma_tile_mmain(A, k);
                        int ldbk = plasma_tile_mmain(B, k);
                        double zbeta = k == 0 ? beta : 1.0;
                        if (k < m) {
                            core_omp_dgemm(
                                PlasmaNoTrans, PlasmaNoTrans,
                                mvcm, nvcn, mvck,
                                alpha, A(m, k), ldam,
                                       B(k, n), ldbk,
                                zbeta, C(m, n), ldcm,
                                sequence, request);
                        }
                        else {
                            if (k == m) {
                                core_omp_dsymm(
                                    side, uplo,
                                    mvcm, nvcn,
                                    alpha, A(k, k), ldak,
                                           B(k, n), ldbk,
                                    zbeta, C(m, n), ldcm,
                                    sequence, request);
                            }
                            else {
                                core_omp_dgemm(
                                    PlasmaTrans, PlasmaNoTrans,
                                    mvcm, nvcn, mvck,
                                    alpha, A(k, m), ldak,
                                           B(k, n), ldbk,
                                    zbeta, C(m, n), ldcm,
                                    sequence, request);
                            }
                        }
                    }
                }
                //===========================
                // PlasmaLeft / PlasmaUpper
                //===========================
                else {
                    for (int k = 0; k < C.mt; k++) {
                        int mvck = plasma_tile_mview(C, k);
                        int ldak = plasma_tile_mmain(A, k);
                        int ldbk = plasma_tile_mmain(B, k);
                        double zbeta = k == 0 ? beta : 1.0;
                        if (k < m) {
                            core_omp_dgemm(
                                PlasmaTrans, PlasmaNoTrans,
                                mvcm, nvcn, mvck,
                                alpha, A(k, m), ldak,
                                       B(k, n), ldbk,
                                zbeta, C(m, n), ldcm,
                                sequence, request);
                        }
                        else {
                            if (k == m) {
                                core_omp_dsymm(
                                    side, uplo,
                                    mvcm, nvcn,
                                    alpha, A(k, k), ldak,
                                           B(k, n), ldbk,
                                    zbeta, C(m, n), ldcm,
                                    sequence, request);
                            }
                            else {
                                core_omp_dgemm(
                                    PlasmaNoTrans, PlasmaNoTrans,
                                    mvcm, nvcn, mvck,
                                    alpha, A(m, k), ldam,
                                           B(k, n), ldbk,
                                    zbeta, C(m, n), ldcm,
                                    sequence, request);
                            }
                        }
                    }
                }
            }
            else {
                int ldan = plasma_tile_mmain(A, n);
                int ldbm = plasma_tile_mmain(B, m);
                //============================
                // PlasmaRight / PlasmaLower
                //============================
                if (uplo == PlasmaLower) {
                    for (int k = 0; k < C.nt; k++) {
                        int nvck = plasma_tile_nview(C, k);
                        int ldak   = plasma_tile_mmain(A, k);
                        double zbeta = k == 0 ? beta : 1.0;
                        if (k < n) {
                            core_omp_dgemm(
                                PlasmaNoTrans, PlasmaTrans,
                                mvcm, nvcn, nvck,
                                alpha, B(m, k), ldbm,
                                       A(n, k), ldan,
                                zbeta, C(m, n), ldcm,
                                sequence, request);
                        }
                        else {
                            if (n == k) {
                                core_omp_dsymm(
                                    side, uplo,
                                    mvcm, nvcn,
                                    alpha, A(k, k), ldak,
                                           B(m, k), ldbm,
                                    zbeta, C(m, n), ldcm,
                                    sequence, request);
                            }
                            else {
                                core_omp_dgemm(
                                    PlasmaNoTrans, PlasmaNoTrans,
                                    mvcm, nvcn, nvck,
                                    alpha, B(m, k), ldbm,
                                           A(k, n), ldak,
                                    zbeta, C(m, n), ldcm,
                                    sequence, request);
                            }
                        }
                    }
                }
                //============================
                // PlasmaRight / PlasmaUpper
                //============================
                else {
                    for (int k = 0; k < C.nt; k++) {
                        int nvck = plasma_tile_nview(C, k);
                        int ldak = plasma_tile_mmain(A, k);
                        double zbeta = k == 0 ? beta : 1.0;
                        if (k < n) {
                            core_omp_dgemm(
                                PlasmaNoTrans, PlasmaNoTrans,
                                mvcm, nvcn, nvck,
                                alpha, B(m, k), ldbm,
                                       A(k, n), ldak,
                                zbeta, C(m, n), ldcm,
                                sequence, request);
                        }
                        else {
                            if (n == k) {
                                core_omp_dsymm(
                                    side, uplo,
                                    mvcm, nvcn,
                                    alpha, A(k, k), ldak,
                                           B(m, k), ldbm,
                                    zbeta, C(m, n), ldcm,
                                    sequence, request);
                            }
                            else {
                                core_omp_dgemm(
                                    PlasmaNoTrans, PlasmaTrans,
                                    mvcm, nvcn, nvck,
                                    alpha, B(m, k), ldbm,
                                           A(n, k), ldan,
                                    zbeta, C(m, n), ldcm,
                                    sequence, request);
                            }
                        }
                    }
                }
            }
        }
    }
}
