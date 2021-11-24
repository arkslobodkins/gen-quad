/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from compute/pclag2z.c, mixed zc -> ds, Mon Nov 22 19:21:35 2021
 *
 **/

#include "plasma_async.h"
#include "plasma_descriptor.h"
#include "plasma_types.h"
#include "plasma_internal.h"
#include "core_blas_ds.h"

#define As(m, n) (float*)plasma_tile_addr(As, m, n)
#define  A(m, n) (double*)plasma_tile_addr( A, m, n)

/***************************************************************************//**
 * Parallel tile conversion of matrix precision from single complex to
 * double complex.
 * @see plasma_omp_slag2d
 ******************************************************************************/
void plasma_pslag2d(plasma_desc_t As, plasma_desc_t A,
                    plasma_sequence_t *sequence, plasma_request_t *request)
{
    // Return if failed sequence.
    if (sequence->status != PlasmaSuccess)
        return;

    for (int m = 0; m < As.mt; m++) {
        int am  = plasma_tile_mview(As, m);
        int lda = plasma_tile_mmain(As, m);
        int ldb = plasma_tile_mmain(A,  m);
        for (int n = 0; n < As.nt; n++) {
            int an = plasma_tile_nview(As, n);
            core_omp_slag2d(
                am, an,
                As(m, n), lda,
                A(m, n),  ldb,
                sequence, request);
        }
    }
}
