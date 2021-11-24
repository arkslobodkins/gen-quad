/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from compute/pzlacpy.c, normal z -> c, Mon Nov 22 19:21:08 2021
 *
 **/

#include "plasma_async.h"
#include "plasma_descriptor.h"
#include "plasma_internal.h"
#include "plasma_types.h"
#include "core_blas.h"

#define A(m, n) (plasma_complex32_t*)plasma_tile_addr(A, m, n)
#define B(m, n) (plasma_complex32_t*)plasma_tile_addr(B, m, n)

/***************************************************************************//**
 * Parallel tile matrix copy.
 * @see plasma_omp_clacpy
 ******************************************************************************/
void plasma_pclacpy(plasma_enum_t uplo, plasma_desc_t A, plasma_desc_t B,
                    plasma_sequence_t *sequence, plasma_request_t *request)
{
    // Return if failed sequence.
    if (sequence->status != PlasmaSuccess)
        return;

    switch (uplo) {
    //==============
    // PlasmaUpper
    //==============
    case PlasmaUpper:
        for (int m = 0; m < A.mt; m++) {
            int mvam = plasma_tile_mview(A, m);
            int ldam = plasma_tile_mmain(A, m);
            int ldbm = plasma_tile_mmain(B, m);
            if (m < A.nt) {
                int nvan = plasma_tile_nview(A, m);
                core_omp_clacpy(
                    PlasmaUpper,
                    mvam, nvan,
                    A(m, m), ldam,
                    B(m, m), ldbm,
                    sequence, request);
            }
            for (int n = m+1; n < A.nt; n++) {
                int nvan = plasma_tile_nview(A, n);
                core_omp_clacpy(
                    PlasmaGeneral,
                    mvam, nvan,
                    A(m, n), ldam,
                    B(m, n), ldbm,
                    sequence, request);
            }
        }
        break;
    //==============
    // PlasmaLower
    //==============
    case PlasmaLower:
        for (int m = 0; m < A.mt; m++) {
            int mvam = plasma_tile_mview(A, m);
            int ldam = plasma_tile_mmain(A, m);
            int ldbm = plasma_tile_mmain(B, m);
            if (m < A.nt) {
                int nvan = plasma_tile_nview(A, m);
                core_omp_clacpy(
                    PlasmaLower,
                    mvam, nvan,
                    A(m, m), ldam,
                    B(m, m), ldbm,
                    sequence, request);
            }
            for (int n = 0; n < imin(m, A.nt); n++) {
                int nvan = plasma_tile_nview(A, n);
                core_omp_clacpy(
                    PlasmaGeneral,
                    mvam, nvan,
                    A(m, n), ldam,
                    B(m, n), ldbm,
                    sequence, request);
            }
        }
        break;
    //================
    // PlasmaGeneral
    //================
    case PlasmaGeneral:
    default:
        for (int m = 0; m < A.mt; m++) {
            int mvam = plasma_tile_mview(A, m);
            int ldam = plasma_tile_mmain(A, m);
            int ldbm = plasma_tile_mmain(B, m);
            for (int n = 0; n < A.nt; n++) {
                int nvan = plasma_tile_nview(A, n);
                core_omp_clacpy(
                    PlasmaGeneral,
                    mvam, nvan,
                    A(m, n), ldam,
                    B(m, n), ldbm,
                    sequence, request);
            }
        }
    }
}
