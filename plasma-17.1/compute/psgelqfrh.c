/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from compute/pzgelqfrh.c, normal z -> s, Thu Mar 10 18:57:37 2022
 *
 **/

#include "plasma_async.h"
#include "plasma_context.h"
#include "plasma_descriptor.h"
#include "plasma_types.h"
#include "plasma_internal.h"
#include "plasma_rh_tree.h"
#include "core_blas_s.h"

#define A(m, n) (float*)plasma_tile_addr(A, m, n)
#define T(m, n) (float*)plasma_tile_addr(T, m, n)
#define T2(m, n) (float*)plasma_tile_addr(T, m, n+(T.nt/2))
/***************************************************************************//**
 *  Parallel tile LQ factorization based on a tree Householder reduction
 * @see plasma_omp_sgelqf
 **/
void plasma_psgelqfrh(plasma_desc_t A, plasma_desc_t T,
                      plasma_workspace_t work,
                      plasma_sequence_t *sequence, plasma_request_t *request)
{
    // Return if failed sequence.
    if (sequence->status != PlasmaSuccess)
        return;

    // Precompute order of LQ operations - compute it as for QR
    // and transpose it.
    int *operations = NULL;
    int num_operations;
    // Transpose m and n to reuse the QR tree.
    plasma_rh_tree_operations(A.nt, A.mt, &operations, &num_operations);

    // Set inner blocking from the T tile row-dimension.
    int ib = T.mb;

    for (int iop = 0; iop < num_operations; iop++) {
        int j, k, kpiv;
        plasma_enum_t kernel;
        // j is row, k and kpiv are columns
        plasma_rh_tree_get_operation(operations, iop, &kernel, &j, &k, &kpiv);

        int mvaj    = plasma_tile_mview(A, j);
        int nvak    = plasma_tile_nview(A, k);
        int ldaj    = plasma_tile_mmain(A, j);

        if (kernel == PlasmaGeKernel) {
            // triangularization
            core_omp_sgelqt(
                mvaj, nvak, ib,
                A(j, k), ldaj,
                T(j, k), T.mb,
                work,
                sequence, request);

            for (int jj = j + 1; jj < A.mt; jj++) {
                int mvajj = plasma_tile_mview(A, jj);
                int ldajj = plasma_tile_mmain(A, jj);

                core_omp_sormlq(
                    PlasmaRight, PlasmaTrans,
                    mvajj, nvak, imin(mvaj, nvak), ib,
                    A(j,  k), ldaj,
                    T(j,  k), T.mb,
                    A(jj, k), ldajj,
                    work,
                    sequence, request);
            }
        }
        else if (kernel == PlasmaTtKernel) {
            // elimination of the tile
            int nvakpiv = plasma_tile_nview(A, kpiv);

            core_omp_sttlqt(
                mvaj, nvak, ib,
                A(j,  kpiv), ldaj,
                A(j,  k),    ldaj,
                T2(j, k),    T.mb,
                work,
                sequence, request);

            for (int jj = j + 1; jj < A.mt; jj++) {
                int mvajj = plasma_tile_mview(A, jj);
                int ldajj = plasma_tile_mmain(A, jj);

                core_omp_sttmlq(
                    PlasmaRight, PlasmaTrans,
                    mvajj, nvakpiv, mvajj, nvak, imin(mvaj, nvakpiv+nvak), ib,
                    A(jj, kpiv), ldajj,
                    A(jj, k),    ldajj,
                    A(j,  k),    ldaj,
                    T2(j, k),    T.mb,
                    work,
                    sequence, request);
            }
        }
        else if (kernel == PlasmaTsKernel) {
            // elimination of the tile
            int nvakpiv = plasma_tile_nview(A, kpiv);

            core_omp_stslqt(
                mvaj, nvak, ib,
                A(j,  kpiv), ldaj,
                A(j,  k),    ldaj,
                T2(j, k),    T.mb,
                work,
                sequence, request);

            for (int jj = j + 1; jj < A.mt; jj++) {
                int mvajj = plasma_tile_mview(A, jj);
                int ldajj = plasma_tile_mmain(A, jj);

                core_omp_stsmlq(
                    PlasmaRight, PlasmaTrans,
                    mvajj, nvakpiv, mvajj, nvak, imin(mvaj, nvakpiv+nvak), ib,
                    A(jj, kpiv), ldajj,
                    A(jj, k),    ldajj,
                    A(j,  k),    ldaj,
                    T2(j, k),    T.mb,
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
