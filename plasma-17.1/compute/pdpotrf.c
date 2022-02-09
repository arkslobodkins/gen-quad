/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from compute/pzpotrf.c, normal z -> d, Tue Feb  8 19:16:04 2022
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

/***************************************************************************//**
 *  Parallel tile Cholesky factorization.
 * @see plasma_omp_dpotrf
 ******************************************************************************/
void plasma_pdpotrf(plasma_enum_t uplo, plasma_desc_t A,
                    plasma_sequence_t *sequence, plasma_request_t *request)
{
    // Return if failed sequence.
    if (sequence->status != PlasmaSuccess)
        return;

    //==============
    // PlasmaLower
    //==============
    if (uplo == PlasmaLower) {
        for (int k = 0; k < A.mt; k++) {
            int mvak = plasma_tile_mview(A, k);
            int ldak = plasma_tile_mmain(A, k);
            core_omp_dpotrf(
                PlasmaLower, mvak,
                A(k, k), ldak,
                A.nb*k,
                sequence, request);

            for (int m = k+1; m < A.mt; m++) {
                int mvam = plasma_tile_mview(A, m);
                int ldam = plasma_tile_mmain(A, m);
                core_omp_dtrsm(
                    PlasmaRight, PlasmaLower,
                    PlasmaConjTrans, PlasmaNonUnit,
                    mvam, A.mb,
                    1.0, A(k, k), ldak,
                         A(m, k), ldam,
                    sequence, request);
            }
            for (int m = k+1; m < A.mt; m++) {
                int mvam = plasma_tile_mview(A, m);
                int ldam = plasma_tile_mmain(A, m);
                core_omp_dsyrk(
                    PlasmaLower, PlasmaNoTrans,
                    mvam, A.mb,
                    -1.0, A(m, k), ldam,
                     1.0, A(m, m), ldam,
                    sequence, request);

                for (int n = k+1; n < m; n++) {
                    int ldan = plasma_tile_mmain(A, n);
                    core_omp_dgemm(
                        PlasmaNoTrans, PlasmaConjTrans,
                        mvam, A.mb, A.mb,
                        -1.0, A(m, k), ldam,
                              A(n, k), ldan,
                         1.0, A(m, n), ldam,
                        sequence, request);
                }
            }
        }
    }
    //==============
    // PlasmaUpper
    //==============
    else {
        for (int k = 0; k < A.nt; k++) {
            int nvak = plasma_tile_nview(A, k);
            int ldak = plasma_tile_mmain(A, k);
            core_omp_dpotrf(
                PlasmaUpper, nvak,
                A(k, k), ldak,
                A.nb*k,
                sequence, request);

            for (int m = k+1; m < A.nt; m++) {
                int nvam = plasma_tile_nview(A, m);
                core_omp_dtrsm(
                    PlasmaLeft, PlasmaUpper,
                    PlasmaConjTrans, PlasmaNonUnit,
                    A.nb, nvam,
                    1.0, A(k, k), ldak,
                         A(k, m), ldak,
                    sequence, request);
            }
            for (int m = k+1; m < A.nt; m++) {
                int nvam = plasma_tile_nview(A, m);
                int ldam = plasma_tile_mmain(A, m);
                core_omp_dsyrk(
                    PlasmaUpper, PlasmaConjTrans,
                    nvam, A.mb,
                    -1.0, A(k, m), ldak,
                     1.0, A(m, m), ldam,
                    sequence, request);

                for (int n = k+1; n < m; n++) {
                    int ldan = plasma_tile_mmain(A, n);
                    core_omp_dgemm(
                        PlasmaConjTrans, PlasmaNoTrans,
                        A.mb, nvam, A.mb,
                        -1.0, A(k, n), ldak,
                              A(k, m), ldak,
                         1.0, A(n, m), ldan,
                        sequence, request);
                }
            }
        }
    }
}
