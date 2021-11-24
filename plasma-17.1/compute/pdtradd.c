/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from compute/pztradd.c, normal z -> d, Mon Nov 22 19:21:40 2021
 *
 **/

#include "plasma_async.h"
#include "plasma_context.h"
#include "plasma_descriptor.h"
#include "plasma_internal.h"
#include "plasma_types.h"
#include "plasma_workspace.h"
#include "core_blas.h"

#define A(m,n) (double*)plasma_tile_addr(A, m, n)
#define B(m,n) (double*)plasma_tile_addr(B, m, n)

/***************************************************************************//**
 * Parallel tile matrix-matrix addition.
 * @see plasma_omp_dtradd
 ******************************************************************************/
void plasma_pdtradd(plasma_enum_t uplo, plasma_enum_t transa,
                    double alpha,  plasma_desc_t A,
                    double beta,   plasma_desc_t B,
                    plasma_sequence_t *sequence, plasma_request_t *request)
{
    // Return if failed sequence.
    if (sequence->status != PlasmaSuccess)
        return;

    switch (uplo) {
    case PlasmaLower:
        //==============================
        // PlasmaLower / PlasmaNoTrans
        //==============================
        if (transa == PlasmaNoTrans) {
            for (int n = 0; n < imin(B.mt,B.nt); n++) {
                int mvbn = plasma_tile_mview(B, n);
                int nvbn = plasma_tile_nview(B, n);
                int ldan = plasma_tile_mmain(A, n);
                int ldbn = plasma_tile_mmain(B, n);
                core_omp_dtradd(
                    uplo, transa,
                    mvbn, nvbn,
                    alpha, A(n, n), ldan,
                    beta,  B(n, n), ldbn,
                    sequence, request);

                for (int m = n+1; m < B.mt; m++) {
                    int mvbm = plasma_tile_mview(B, m);
                    int ldam = plasma_tile_mmain(A, m);
                    int ldbm = plasma_tile_mmain(B, m);
                    core_omp_dgeadd(
                        transa,
                        mvbm, nvbn,
                        alpha, A(m, n), ldam,
                        beta,  B(m, n), ldbm,
                        sequence, request);
                }
            }
        }
        //===================================
        // PlasmaLower / Plasma[_Conj]Trans
        //===================================
        else {
            for (int n = 0; n < imin(B.mt,B.nt); n++) {
                int mvbn = plasma_tile_mview(B, n);
                int nvbn = plasma_tile_nview(B, n);
                int ldan = plasma_tile_mmain(A, n);
                int ldbn = plasma_tile_mmain(B, n);
                core_omp_dtradd(
                    uplo, transa,
                    mvbn, nvbn,
                    alpha, A(n, n), ldan,
                    beta,  B(n, n), ldbn,
                    sequence, request);

                for (int m = n+1; m < B.mt; m++) {
                    int mvbm = plasma_tile_mview(B, m);
                    int ldbm = plasma_tile_mmain(B, m);
                    core_omp_dgeadd(
                        transa,
                        mvbm, nvbn,
                        alpha, A(n, m), ldan,
                        beta,  B(m, n), ldbm,
                        sequence, request);
                }
            }
        }
        break;
    case PlasmaUpper:
        //==============================
        // PlasmaUpper / PlasmaNoTrans
        //==============================
        if (transa == PlasmaNoTrans) {
            for (int m = 0; m < imin(B.mt,B.nt); m++) {
                int mvbm = plasma_tile_mview(B, m);
                int nvbm = plasma_tile_nview(B, m);
                int ldam = plasma_tile_mmain(A, m);
                int ldbm = plasma_tile_mmain(B, m);
                core_omp_dtradd(
                    uplo, transa,
                    mvbm, nvbm,
                    alpha, A(m, m), ldam,
                    beta,  B(m, m), ldbm,
                    sequence, request);

                for (int n = m+1; n < B.nt; n++) {
                    int nvbn = plasma_tile_nview(B, n);
                    core_omp_dgeadd(
                        transa,
                        mvbm, nvbn,
                        alpha, A(m, n), ldam,
                        beta,  B(m, n), ldbm,
                        sequence, request);
                }
            }
        }
        //===================================
        // PlasmaUpper / Plasma[_Conj]Trans
        //===================================
        else {
            for (int m = 0; m < imin(B.mt,B.nt); m++) {
                int mvbm = plasma_tile_mview(B, m);
                int nvbm = plasma_tile_nview(B, m);
                int ldam = plasma_tile_mmain(A, m);
                int ldbm = plasma_tile_mmain(B, m);
                core_omp_dtradd(
                    uplo, transa,
                    mvbm, nvbm,
                    alpha, A(m, m), ldam,
                    beta,  B(m, m), ldbm,
                    sequence, request);

                for (int n = m+1; n < B.nt; n++) {
                    int nvbn = plasma_tile_nview(B, n);
                    int ldan = plasma_tile_mmain(A, n);
                    core_omp_dgeadd(
                        transa,
                        mvbm, nvbn,
                        alpha, A(n, m), ldan,
                        beta,  B(m, n), ldbm,
                        sequence, request);
                }
            }
        }
        break;
    }
}
