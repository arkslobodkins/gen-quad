/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from compute/pdzamax.c, normal z -> s, Mon Nov 22 19:21:42 2021
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

/******************************************************************************/
void plasma_psamax(plasma_enum_t colrow,
                    plasma_desc_t A, float *work, float *values,
                    plasma_sequence_t *sequence, plasma_request_t *request)
{
    // Return if failed sequence.
    if (sequence->status != PlasmaSuccess)
        return;

    switch (colrow) {
    //===================
    // PlasmaColumnwise
    //===================
    case PlasmaColumnwise:
        for (int m = 0; m < A.mt; m++) {
            int mvam = plasma_tile_mview(A, m);
            int ldam = plasma_tile_mmain(A, m);
            for (int n = 0; n < A.nt; n++) {
                int nvan = plasma_tile_nview(A, n);
                core_omp_samax(PlasmaColumnwise,
                                mvam, nvan,
                                A(m, n), ldam,
                                &work[A.n*m+n*A.nb],
                                sequence, request);
            }
        }
        #pragma omp taskwait
        core_omp_samax(PlasmaRowwise,
                       A.n, A.mt,
                       work, A.n,
                       values,
                       sequence, request);
        break;
    //================
    // PlasmaRowwise
    //================
    case PlasmaRowwise:
        for (int m = 0; m < A.mt; m++) {
            int mvam = plasma_tile_mview(A, m);
            int ldam = plasma_tile_mmain(A, m);
            for (int n = 0; n < A.nt; n++) {
                int nvan = plasma_tile_nview(A, n);
                core_omp_samax(PlasmaRowwise,
                                mvam, nvan,
                                A(m, n), ldam,
                                &work[A.m*n+m*A.mb],
                                sequence, request);
            }
        }
        #pragma omp taskwait
        core_omp_samax(PlasmaRowwise,
                       A.m, A.nt,
                       work, A.m,
                       values,
                       sequence, request);
    }
}
