/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from compute/pzgeqrfrh.c, normal z -> d, Thu Mar 10 18:57:46 2022
 *
 **/

#include "plasma_async.h"
#include "plasma_context.h"
#include "plasma_descriptor.h"
#include "plasma_types.h"
#include "plasma_internal.h"
#include "plasma_rh_tree.h"
#include "core_blas_d.h"

#define A(m, n) (double*)plasma_tile_addr(A, m, n)
#define T(m, n) (double*)plasma_tile_addr(T, m, n)
#define T2(m, n) (double*)plasma_tile_addr(T, m, n+(T.nt/2))
/***************************************************************************//**
 *  Parallel tile QR factorization based on a tree Householder reduction
 * @see plasma_omp_dgeqrf
 **/
void plasma_pdgeqrfrh(plasma_desc_t A, plasma_desc_t T,
                      plasma_workspace_t work,
                      plasma_sequence_t *sequence, plasma_request_t *request)
{
    // Return if failed sequence.
    if (sequence->status != PlasmaSuccess)
        return;

    // Precompute order of QR operations.
    int *operations = NULL;
    int num_operations;
    plasma_rh_tree_operations(A.mt, A.nt, &operations, &num_operations);

    // Set inner blocking from the T tile row-dimension.
    int ib = T.mb;

    for (int iop = 0; iop < num_operations; iop++) {
        int j, k, kpiv;
        plasma_enum_t kernel;
        plasma_rh_tree_get_operation(operations, iop, &kernel, &j, &k, &kpiv);

        int nvaj    = plasma_tile_nview(A, j);
        int mvak    = plasma_tile_mview(A, k);
        int ldak    = plasma_tile_mmain(A, k);

        if (kernel == PlasmaGeKernel) {
            // triangularization
            core_omp_dgeqrt(
                mvak, nvaj, ib,
                A(k, j), ldak,
                T(k, j), T.mb,
                work,
                sequence, request);

            for (int jj = j + 1; jj < A.nt; jj++) {
                int nvajj = plasma_tile_nview(A, jj);

                core_omp_dormqr(
                    PlasmaLeft, PlasmaTrans,
                    mvak, nvajj, imin(mvak, nvaj), ib,
                    A(k, j), ldak,
                    T(k, j), T.mb,
                    A(k, jj), ldak,
                    work,
                    sequence, request);
            }
        }
        else if (kernel == PlasmaTtKernel) {
            // elimination of the tile
            int mvakpiv = plasma_tile_mview(A, kpiv);
            int ldakpiv = plasma_tile_mmain(A, kpiv);

            core_omp_dttqrt(
                mvak, nvaj, ib,
                A(kpiv, j), ldakpiv,
                A(k,  j),   ldak,
                T2(k, j),   T.mb,
                work,
                sequence, request);

            for (int jj = j + 1; jj < A.nt; jj++) {
                int nvajj = plasma_tile_nview(A, jj);

                core_omp_dttmqr(
                    PlasmaLeft, PlasmaTrans,
                    mvakpiv, nvajj, mvak, nvajj, imin(mvakpiv+mvak, nvaj), ib,
                    A(kpiv, jj), ldakpiv,
                    A(k,    jj), ldak,
                    A(k,    j),  ldak,
                    T2(k,   j),  T.mb,
                    work,
                    sequence, request);
            }
        }
        else if (kernel == PlasmaTsKernel) {
            // elimination of the tile
            int mvakpiv = plasma_tile_mview(A, kpiv);
            int ldakpiv = plasma_tile_mmain(A, kpiv);

            core_omp_dtsqrt(
                mvak, nvaj, ib,
                A(kpiv, j), ldakpiv,
                A(k,  j),   ldak,
                T2(k, j),   T.mb,
                work,
                sequence, request);

            for (int jj = j + 1; jj < A.nt; jj++) {
                int nvajj = plasma_tile_nview(A, jj);

                core_omp_dtsmqr(
                    PlasmaLeft, PlasmaTrans,
                    mvakpiv, nvajj, mvak, nvajj, imin(mvakpiv+mvak, nvaj), ib,
                    A(kpiv, jj), ldakpiv,
                    A(k,    jj), ldak,
                    A(k,    j),  ldak,
                    T2(k,   j),  T.mb,
                    work,
                    sequence, request);
            }
        }
        else {
            plasma_error("illegal kernel");
            plasma_request_fail(sequence, request, PlasmaErrorIllegalValue);
        }
    }

    free(operations);
}
