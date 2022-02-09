/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from compute/pzlansy.c, normal z -> s, Tue Feb  8 19:15:52 2022
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

/***************************************************************************//**
 *  Parallel tile calculation of max, one, infinity or Frobenius matrix norm
 *  for a symmetric matrix.
 ******************************************************************************/
void plasma_pslansy(plasma_enum_t norm, plasma_enum_t uplo,
                    plasma_desc_t A, float *work, float *value,
                    plasma_sequence_t *sequence, plasma_request_t *request)
{
    // Return if failed sequence.
    if (sequence->status != PlasmaSuccess)
        return;

    switch (norm) {
    float stub;
    float *workspace;
    float *scale;
    float *sumsq;
    //================
    // PlasmaMaxNorm
    //================
    case PlasmaMaxNorm:
        for (int m = 0; m < A.mt; m++) {
            int mvam = plasma_tile_mview(A, m);
            int ldam = plasma_tile_mmain(A, m);
            if (uplo == PlasmaLower) {
                for (int n = 0; n < m; n++) {
                    int nvan = plasma_tile_nview(A, n);
                    core_omp_slange(PlasmaMaxNorm,
                                    mvam, nvan,
                                    A(m, n), ldam,
                                    &stub, &work[A.mt*n+m],
                                    sequence, request);
                }
            }
            else { // PlasmaUpper
                for (int n = m+1; n < A.nt; n++) {
                    int nvan = plasma_tile_nview(A, n);
                    core_omp_slange(PlasmaMaxNorm,
                                    mvam, nvan,
                                    A(m, n), ldam,
                                    &stub, &work[A.mt*n+m],
                                    sequence, request);
                }
            }
            core_omp_slansy(PlasmaMaxNorm, uplo,
                            mvam,
                            A(m, m), ldam,
                            &stub, &work[A.mt*m+m],
                            sequence, request);
        }
        #pragma omp taskwait
        core_omp_slansy(PlasmaMaxNorm, uplo,
                        A.nt,
                        work, A.mt,
                        &stub, value,
                        sequence, request);
        break;
    //================
    // PlasmaOneNorm
    //================
    case PlasmaOneNorm:
    case PlasmaInfNorm:
        for (int m = 0; m < A.mt; m++) {
            int mvam = plasma_tile_mview(A, m);
            int ldam = plasma_tile_mmain(A, m);
            if (uplo == PlasmaLower) {
                for (int n = 0; n < m; n++) {
                    int nvan = plasma_tile_nview(A, n);
                    core_omp_slange_aux(PlasmaOneNorm,
                                        mvam, nvan,
                                        A(m, n), ldam,
                                        &work[A.n*m+n*A.nb],
                                        sequence, request);
                    core_omp_slange_aux(PlasmaInfNorm,
                                        mvam, nvan,
                                        A(m, n), ldam,
                                        &work[A.n*n+m*A.nb],
                                        sequence, request);
                }
            }
            else { // PlasmaUpper
                for (int n = m+1; n < A.nt; n++) {
                    int nvan = plasma_tile_nview(A, n);
                    core_omp_slange_aux(PlasmaOneNorm,
                                        mvam, nvan,
                                        A(m, n), ldam,
                                        &work[A.n*m+n*A.nb],
                                        sequence, request);
                    core_omp_slange_aux(PlasmaInfNorm,
                                        mvam, nvan,
                                        A(m, n), ldam,
                                        &work[A.n*n+m*A.nb],
                                        sequence, request);
                }
            }
            core_omp_slansy_aux(PlasmaOneNorm, uplo,
                                mvam,
                                A(m, m), ldam,
                                &work[A.n*m+m*A.nb],
                                sequence, request);
        }
        #pragma omp taskwait
        workspace = work + A.mt*A.n;
        core_omp_slange(PlasmaInfNorm,
                        A.n, A.mt,
                        work, A.n,
                        workspace, value,
                        sequence, request);
        break;
    //======================
    // PlasmaFrobeniusNorm
    //======================
    case PlasmaFrobeniusNorm:
        scale = work;
        sumsq = work + A.mt*A.nt;
        for (int m = 0; m < A.mt; m++) {
            int mvam = plasma_tile_mview(A, m);
            int ldam = plasma_tile_mmain(A, m);
            if (uplo == PlasmaLower) {
                for (int n = 0; n < m; n++) {
                    int nvan = plasma_tile_nview(A, n);
                    core_omp_sgessq(mvam, nvan,
                                    A(m, n), ldam,
                                    &scale[A.mt*n+m], &sumsq[A.mt*n+m],
                                    sequence, request);
                }
            }
            else { // PlasmaUpper
                for (int n = m+1; n < A.nt; n++) {
                    int nvan = plasma_tile_nview(A, n);
                    core_omp_sgessq(mvam, nvan,
                                    A(m, n), ldam,
                                    &scale[A.mt*m+n], &sumsq[A.mt*m+n],
                                    sequence, request);
                }
            }
            core_omp_ssyssq(uplo,
                            mvam,
                            A(m, m), ldam,
                            &scale[A.mt*m+m], &sumsq[A.mt*m+m],
                            sequence, request);
        }
        #pragma omp taskwait
        core_omp_ssyssq_aux(A.mt, A.nt,
                            scale, sumsq,
                            value,
                            sequence, request);
        break;
    }
}
