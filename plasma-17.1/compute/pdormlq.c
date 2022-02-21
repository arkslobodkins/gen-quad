/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from compute/pzunmlq.c, normal z -> d, Sun Feb 20 23:41:19 2022
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
#define T(m, n) (double*)plasma_tile_addr(T, m, n)
#define B(m, n) (double*)plasma_tile_addr(B, m, n)

/***************************************************************************//**
 *  Parallel application of Q using tile V - LQ factorization
 * @see plasma_omp_dgelqs
 **/
void plasma_pdormlq(plasma_enum_t side, plasma_enum_t trans,
                    plasma_desc_t A, plasma_desc_t T, plasma_desc_t B,
                    plasma_workspace_t work,
                    plasma_sequence_t *sequence, plasma_request_t *request)
{
    // Return if failed sequence.
    if (sequence->status != PlasmaSuccess)
        return;

    // Set inner blocking from the T tile row-dimension.
    int ib = T.mb;

    if (side == PlasmaLeft) {
        //=============================
        // PlasmaLeft / PlasmaNoTrans
        //=============================
        if (trans == PlasmaNoTrans) {
            for (int k = 0; k < imin(A.mt, A.nt); k++) {
                int mvbk = plasma_tile_mview(B, k);
                int mvak = plasma_tile_mview(A, k);
                int nvak = plasma_tile_nview(A, k);
                int ldak = plasma_tile_mmain(A, k);
                int ldbk = plasma_tile_mmain(B, k);
                for (int n = 0; n < B.nt; n++) {
                    int nvbn = plasma_tile_nview(B, n);
                    core_omp_dormlq(
                            side, trans,
                            mvbk, nvbn, imin(nvak, mvak), ib,
                            A(k, k), ldak,
                            T(k, k), T.mb,
                            B(k, n), ldbk,
                            work,
                            sequence, request);
                }
                for (int m = k+1; m < B.mt; m++) {
                    int mvbm = plasma_tile_mview(B, m);
                    int ldbm = plasma_tile_mmain(B, m);
                    for (int n = 0; n < B.nt; n++) {
                        int nvbn = plasma_tile_nview(B, n);
                        core_omp_dtsmlq(
                                side, trans,
                                B.mb, nvbn, mvbm, nvbn, mvak, ib,
                                B(k, n), ldbk,
                                B(m, n), ldbm,
                                A(k, m), ldak,
                                T(k, m), T.mb,
                                work,
                                sequence, request);
                    }
                }
            }
        }
        //==================================
        // PlasmaLeft / Plasma[_Conj]Trans
        //==================================
        else {
            for (int k = imin(A.mt, A.nt)-1; k >= 0; k--) {
                int mvbk = plasma_tile_mview(B, k);
                int mvak = plasma_tile_mview(A, k);
                int nvak = plasma_tile_nview(A, k);
                int ldak = plasma_tile_mmain(A, k);
                int ldbk = plasma_tile_mmain(B, k);
                for (int m = B.mt-1; m > k; m--) {
                    int mvbm = plasma_tile_mview(B, m);
                    int ldbm = plasma_tile_mmain(B, m);
                    for (int n = 0; n < B.nt; n++) {
                        int nvbn = plasma_tile_nview(B, n);
                        core_omp_dtsmlq(
                                side, trans,
                                B.mb, nvbn, mvbm, nvbn, mvak, ib,
                                B(k, n), ldbk,
                                B(m, n), ldbm,
                                A(k, m), ldak,
                                T(k, m), T.mb,
                                work,
                                sequence, request);
                    }
                }
                for (int n = 0; n < B.nt; n++) {
                    int nvbn = plasma_tile_nview(B, n);
                    core_omp_dormlq(
                            side, trans,
                            mvbk, nvbn, imin(nvak, mvak), ib,
                            A(k, k), ldak,
                            T(k, k), T.mb,
                            B(k, n), ldbk,
                            work,
                            sequence, request);
                }
            }
        }
    }
    else {
        //==============================
        // PlasmaRight / PlasmaNoTrans
        //==============================
        if (trans == PlasmaNoTrans) {
            for (int k = imin(A.mt, A.nt)-1; k >= 0; k--) {
                int nvbk = plasma_tile_nview(B, k);
                int mvak = plasma_tile_mview(A, k);
                int nvak = plasma_tile_nview(A, k);
                int ldak = plasma_tile_mmain(A, k);
                for (int n = B.nt-1; n > k; n--) {
                    int nvbn = plasma_tile_nview(B, n);
                    for (int m = 0; m < B.mt; m++) {
                        int mvbm = plasma_tile_mview(B, m);
                        int ldbm = plasma_tile_mmain(B, m);
                        core_omp_dtsmlq(
                                side, trans,
                                mvbm, B.nb, mvbm, nvbn, mvak, ib,
                                B(m, k), ldbm,
                                B(m, n), ldbm,
                                A(k, n), ldak,
                                T(k, n), T.mb,
                                work,
                                sequence, request);
                    }
                }
                for (int m = 0; m < B.mt; m++) {
                    int mvbm = plasma_tile_mview(B, m);
                    int ldbm = plasma_tile_mmain(B, m);
                    core_omp_dormlq(
                            side, trans,
                            mvbm, nvbk, imin(nvak, mvak), ib,
                            A(k, k), ldak,
                            T(k, k), T.mb,
                            B(m, k), ldbm,
                            work,
                            sequence, request);
                }
            }
        }
        //===================================
        // PlasmaRight / Plasma[_Conj]Trans
        //===================================
        else {
            for (int k = 0; k < imin(A.mt, A.nt); k++) {
                int nvbk = plasma_tile_nview(B, k);
                int mvak = plasma_tile_mview(A, k);
                int nvak = plasma_tile_nview(A, k);
                int ldak = plasma_tile_mmain(A, k);
                for (int m = 0; m < B.mt; m++) {
                    int mvbm = plasma_tile_mview(B, m);
                    int ldbm = plasma_tile_mmain(B, m);
                    core_omp_dormlq(
                            side, trans,
                            mvbm, nvbk, imin(nvak, mvak), ib,
                            A(k, k), ldak,
                            T(k, k), T.mb,
                            B(m, k), ldbm,
                            work,
                            sequence, request);
                }
                for (int n = k+1; n < B.nt; n++) {
                    int nvbn = plasma_tile_nview(B, n);
                    for (int m = 0; m < B.mt; m++) {
                        int mvbm = plasma_tile_mview(B, m);
                        int ldbm = plasma_tile_mmain(B, m);
                        core_omp_dtsmlq(
                                side, trans,
                                mvbm, B.nb, mvbm, nvbn, mvak, ib,
                                B(m, k), ldbm,
                                B(m, n), ldbm,
                                A(k, n), ldak,
                                T(k, n), T.mb,
                                work,
                                sequence, request);
                    }
                }
            }
        }
    }
}
