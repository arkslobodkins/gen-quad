/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from compute/pzsyr2k.c, normal z -> s, Thu Mar 10 18:58:00 2022
 *
 **/

#include "plasma_async.h"
#include "plasma_context.h"
#include "plasma_descriptor.h"
#include "plasma_internal.h"
#include "plasma_types.h"
#include "plasma_workspace.h"
#include "core_blas.h"

#define A(m, n) (float*)plasma_tile_addr(A, m, n)
#define B(m, n) (float*)plasma_tile_addr(B, m, n)
#define C(m, n) (float*)plasma_tile_addr(C, m, n)

/***************************************************************************//**
 * Parallel tile symmetric rank 2k update.
 * @see plasma_omp_ssyr2k
 ******************************************************************************/
void plasma_pssyr2k(plasma_enum_t uplo, plasma_enum_t trans,
                    float alpha, plasma_desc_t A,
                                              plasma_desc_t B,
                    float beta,  plasma_desc_t C,
                    plasma_sequence_t *sequence, plasma_request_t *request)
{
    // Return if failed sequence.
    if (sequence->status != PlasmaSuccess)
        return;

    for (int n = 0; n < C.nt; n++) {
        int nvcn = plasma_tile_nview(C, n);
        int ldan = plasma_tile_mmain(A, n);
        int ldbn = plasma_tile_mmain(B, n);
        int ldcn = plasma_tile_mmain(C, n);
        //================
        // PlasmaNoTrans
        //================
        if (trans == PlasmaNoTrans) {
            for (int k = 0; k < A.nt; k++) {
                int nvak = plasma_tile_nview(A, k);
                float zbeta = k == 0 ? beta : 1.0;
                core_omp_ssyr2k(
                    uplo, trans,
                    nvcn, nvak,
                    alpha, A(n, k), ldan,
                           B(n, k), ldbn,
                    zbeta, C(n, n), ldcn,
                    sequence, request);
            }
            //==============================
            // PlasmaNoTrans / PlasmaLower
            //==============================
            if (uplo == PlasmaLower) {
                for (int m = n+1; m < C.mt; m++) {
                    int mvcm = plasma_tile_mview(C, m);
                    int ldam = plasma_tile_mmain(A, m);
                    int ldcm = plasma_tile_mmain(C, m);
                    for (int k = 0; k < A.nt; k++) {
                        int nvak = plasma_tile_nview(A, k);
                        float zbeta = k == 0 ? beta : 1.0;
                        core_omp_sgemm(
                            trans, PlasmaTrans,
                            mvcm, nvcn, nvak,
                            alpha, A(m, k), ldam,
                                   B(n, k), ldbn,
                            zbeta, C(m, n), ldcm,
                            sequence, request);
                        core_omp_sgemm(
                            trans, PlasmaTrans,
                            mvcm, nvcn, nvak,
                            alpha, B(m, k), ldam,
                                   A(n, k), ldan,
                            1.0,   C(m, n), ldcm,
                            sequence, request);
                    }
                }
            }
            //==============================
            // PlasmaNoTrans / PlasmaUpper
            //==============================
            else {
                for (int m = n+1; m < C.mt; m++) {
                    int mvcm = plasma_tile_mview(C, m);
                    int ldam = plasma_tile_mmain(A, m);
                    int ldbm = plasma_tile_mmain(B, m);
                    for (int k = 0; k < A.nt; k++) {
                        int nvak = plasma_tile_nview(A, k);
                        float zbeta = k == 0 ? beta : 1.0;
                        core_omp_sgemm(
                            trans, PlasmaTrans,
                            nvcn, mvcm, nvak,
                            alpha, A(n, k), ldan,
                                   B(m, k), ldbm,
                            zbeta, C(n, m), ldcn,
                            sequence, request);
                        core_omp_sgemm(
                            trans, PlasmaTrans,
                            nvcn, mvcm, nvak,
                            alpha, B(n, k), ldan,
                                   A(m, k), ldam,
                            1.0,   C(n, m), ldcn,
                            sequence, request);
                    }
                }
            }
        }
        //=====================
        // Plasma[_Conj]Trans
        //=====================
        else {
            for (int k = 0; k < A.mt; k++) {
                int mvak = plasma_tile_mview(A, k);
                int ldak = plasma_tile_mmain(A, k);
                int ldbk = plasma_tile_mmain(B, k);
                float zbeta = k == 0 ? beta : 1.0;
                core_omp_ssyr2k(
                    uplo, trans,
                    nvcn, mvak,
                    alpha, A(k, n), ldak,
                           B(k, n), ldbk,
                    zbeta, C(n, n), ldcn,
                    sequence, request);
            }
            //===================================
            // Plasma[_Conj]Trans / PlasmaLower
            //===================================
            if (uplo == PlasmaLower) {
                for (int m = n+1; m < C.mt; m++) {
                    int mvcm = plasma_tile_mview(C, m);
                    int ldcm = plasma_tile_mmain(C, m);
                    for (int k = 0; k < A.mt; k++) {
                        int mvak = plasma_tile_mview(A, k);
                        int ldak = plasma_tile_mmain(A, k);
                        int ldbk = plasma_tile_mmain(B, k);
                        float zbeta = k == 0 ? beta : 1.0;
                        core_omp_sgemm(
                            trans, PlasmaNoTrans,
                            mvcm, nvcn, mvak,
                            alpha, A(k, m), ldak,
                                   B(k, n), ldbk,
                            zbeta, C(m, n), ldcm,
                            sequence, request);
                        core_omp_sgemm(
                            trans, PlasmaNoTrans,
                            mvcm, nvcn, mvak,
                            alpha, B(k, m), ldbk,
                                   A(k, n), ldak,
                            1.0,   C(m, n), ldcm,
                            sequence, request);
                    }
                }
            }
            //===================================
            // Plasma[_Conj]Trans / PlasmaUpper
            //===================================
            else {
                for (int m = n+1; m < C.mt; m++) {
                    int mvcm = plasma_tile_mview(C, m);
                    for (int k = 0; k < A.mt; k++) {
                        int mvak = plasma_tile_mview(A, k);
                        int ldak = plasma_tile_mmain(A, k);
                        int ldbk = plasma_tile_mmain(B, k);
                        float zbeta = k == 0 ? beta : 1.0;
                        core_omp_sgemm(
                            trans, PlasmaNoTrans,
                            nvcn, mvcm, mvak,
                            alpha, A(k, n), ldak,
                                   B(k, m), ldbk,
                            zbeta, C(n, m), ldcn,
                            sequence, request);
                        core_omp_sgemm(
                            trans, PlasmaNoTrans,
                            nvcn, mvcm, mvak,
                            alpha, B(k, n), ldbk,
                                   A(k, m), ldak,
                            1.0,   C(n, m), ldcn,
                            sequence, request);
                    }
                }
            }
        }
    }
}
