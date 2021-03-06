/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from compute/pzgbtrf.c, normal z -> d, Thu Mar 10 18:57:41 2022
 *
 **/

#include "plasma_async.h"
#include "plasma_context.h"
#include "plasma_descriptor.h"
#include "plasma_internal.h"
#include "plasma_types.h"
#include "plasma_workspace.h"
#include "core_blas.h"

#define A(m, n) ((double*)plasma_tile_addr(A, m, n))

/******************************************************************************/
void plasma_pdgbtrf(plasma_desc_t A, int *ipiv,
                    plasma_sequence_t *sequence, plasma_request_t *request)
{
    // Return if failed sequence.
    if (sequence->status != PlasmaSuccess)
        return;

    // Read parameters from the context.
    plasma_context_t *plasma = plasma_context_self();
    int ib = plasma->ib;
    int num_panel_threads = plasma->num_panel_threads;
    plasma_barrier_t *barrier = &plasma->barrier;

    for (int k = 0; k < imin(A.mt, A.nt); k++) {
        // for band matrix, gm is a multiple of mb,
        // and there is no a10 submatrix

        int mvak = plasma_tile_mview(A, k);
        int nvak = plasma_tile_nview(A, k);
        int ldak = plasma_tile_mmain_band(A, k, k);

        // panel
        int *ipivk = &ipiv[k*A.mb];
        double *a00 = A(k, k);
        int mak      = imin(A.m-k*A.mb, mvak+A.kl);
        int size_a00 = (A.gm-k*A.mb) * plasma_tile_nmain(A, k);
        int size_i   = imin(mvak, nvak);
        #pragma omp task depend(inout:a00[0:size_a00]) \
                         depend(out:ipivk[0:size_i]) /*\
                         priority(1) */
        {
            if (sequence->status == PlasmaSuccess) {
                for (int rank = 0; rank < num_panel_threads; rank++) {
                    #pragma omp task // priority(1)
                    {
                        // create a view for panel as a "general" submatrix
                        plasma_desc_t view = plasma_desc_view(
                            A, (A.kut-1)*A.mb, k*A.nb, mak, nvak);
                        view.type = PlasmaGeneral;

                        int info = core_dgetrf(view, &ipiv[k*A.mb], ib, rank,
                                               num_panel_threads, barrier);
                        if (info != 0)
                            plasma_request_fail(sequence, request, k*A.mb+info);
                    }
                }
            }
            #pragma omp taskwait
        }
        // update
        // TODO: fills are not tracked, see the one in fork
        for (int n = k+1; n < imin(A.nt, k+A.kut); n++) {
            double *a01 = A(k, n);
            double *a11 = A(k+1, n);

            int nvan = plasma_tile_nview(A, n);
            int size_a01 = ldak*nvan;
            int size_a11 = (A.gm-(k+1)*A.mb)*nvan;

            #pragma omp task depend(in:a00[0:size_a00]) \
                             depend(inout:ipivk[0:size_i]) \
                             depend(inout:a01[0:size_a01]) \
                             depend(inout:a11[0:size_a11]) /*\
                             priority(n == k+1) */
            {
                if (sequence->status == PlasmaSuccess) {
                    // laswp
                    int k1 = k*A.mb+1;
                    int k2 = imin(k*A.mb+A.mb, A.m);
                    plasma_desc_t view =
                        plasma_desc_view(A,
                                        (A.kut-1 + k-n)*A.mb, n*A.nb,
                                         mak, nvan);
                    view.type = PlasmaGeneral;
                    core_dlaswp(
                        PlasmaRowwise, view, 1, k2-k1+1, &ipiv[k*A.mb], 1);

                    // trsm
                    core_dtrsm(PlasmaLeft, PlasmaLower,
                               PlasmaNoTrans, PlasmaUnit,
                               mvak, nvan,
                               1.0, A(k, k), ldak,
                                    A(k, n), plasma_tile_mmain_band(A, k, n));

                    // gemm
                    for (int m = imax(k+1,n-A.kut); m < imin(k+A.klt, A.mt); m++) {
                        int mvam = plasma_tile_mview(A, m);

                        #pragma omp task // priority(n == k+1)
                        {
                            core_dgemm(
                                PlasmaNoTrans, PlasmaNoTrans,
                                mvam, nvan, A.nb,
                                -1.0, A(m, k), plasma_tile_mmain_band(A, m, k),
                                      A(k, n), plasma_tile_mmain_band(A, k, n),
                                1.0,  A(m, n), plasma_tile_mmain_band(A, m, n));
                        }
                    }
                    #pragma omp taskwait
                }
            }
        }
        #pragma omp task depend(in:ipivk[0:size_i])
        if (sequence->status == PlasmaSuccess) {
            if (k > 0) {
                for (int i = 0; i < imin(mak, nvak); i++) {
                    ipiv[k*A.mb+i] += k*A.mb;
                }
            }
        }
    }
}
