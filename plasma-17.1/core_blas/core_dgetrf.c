/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from core_blas/core_zgetrf.c, normal z -> d, Mon Nov 22 19:22:33 2021
 *
 **/

#include "core_blas.h"
#include "core_lapack.h"
#include "plasma_barrier.h"
#include "plasma_descriptor.h"
#include "plasma_internal.h"
#include "plasma_types.h"

#include <omp.h>
#include <assert.h>
#include <math.h>

static int *max_idx;
static double *max_val;
static volatile int info;

#define A(m, n) (double*)plasma_tile_addr(A, m, n)

/******************************************************************************/
int core_dgetrf(plasma_desc_t A, int *ipiv, int ib, int rank, int size,
                plasma_barrier_t *barrier)
{
    // Allocate arrays for parallel max search.
    if (rank == 0) {
        max_idx = (int*)malloc(size*sizeof(int));
        max_val = (double*)malloc(size*sizeof(double));
        info = 0;
    }
    plasma_barrier_wait(barrier);
    assert(max_idx != NULL);
    assert(max_val != NULL);

    double sfmin = LAPACKE_dlamch_work('S');

    for (int k = 0; k < imin(A.m, A.n); k += ib) {
        int kb = imin(imin(A.m, A.n)-k, ib);

        double *a0 = A(0, 0);
        int lda0 = plasma_tile_mmain(A, 0);
        int mva0 = plasma_tile_mview(A, 0);
        int nva0 = plasma_tile_nview(A, 0);

        //======================
        // panel factorization
        //======================
        for (int j = k; j < k+kb; j++) {
            // pivot search
            max_idx[rank] = 0;
            max_val[rank] = a0[j+j*lda0];

            for (int l = rank; l < A.mt; l += size) {
                double *al = A(l, 0);
                int ldal = plasma_tile_mmain(A, l);
                int mval = plasma_tile_mview(A, l);

                if (l == 0) {
                    for (int i = 1; i < mva0-j; i++)
                        if (fabs(a0[j+i+j*lda0]) >
                            fabs(max_val[rank])) {

                            max_val[rank] = a0[j+i+j*lda0];
                            max_idx[rank] = i;
                        }
                }
                else {
                    for (int i = 0; i < mval; i++)
                        if (fabs(al[i+j*ldal]) >
                            fabs(max_val[rank])) {

                            max_val[rank] = al[i+j*ldal];
                            max_idx[rank] = A.mb*l+i-j;
                        }
                }
            }

            plasma_barrier_wait(barrier);
            if (rank == 0)
            {
                // max reduction
                for (int i = 1; i < size; i++) {
                    if (fabs(max_val[i]) >
                        fabs(max_val[0])) {
                        max_val[0] = max_val[i];
                        max_idx[0] = max_idx[i];
                    }
                }

                // pivot adjustment
                int jp = j+max_idx[0];
                ipiv[j] = jp-k+1;

                // singularity check
                if (info == 0 && max_val[0] == 0.0) {
                    info = j+1;
                }
                else {
                    // pivot swap
                    if (jp != j) {
                        double *ap = A(jp/A.mb, 0);
                        int ldap = plasma_tile_mmain(A, jp/A.mb);

                        cblas_dswap(kb,
                                    &a0[j+k*lda0], lda0,
                                    &ap[jp%A.mb+k*ldap], ldap);
                    }
                }
            }
            plasma_barrier_wait(barrier);

            // column scaling and trailing update (all ranks)
            for (int l = rank; l < A.mt; l += size) {
                double *al = A(l, 0);
                int ldal = plasma_tile_mmain(A, l);
                int mval = plasma_tile_mview(A, l);

                if (info == 0) {
                    // column scaling
                    if (fabs(a0[j+j*lda0]) >= sfmin) {
                        if (l == 0) {
                            for (int i = 1; i < mva0-j; i++)
                                a0[j+i+j*lda0] /= a0[j+j*lda0];
                        }
                        else {
                            for (int i = 0; i < mval; i++)
                                al[i+j*ldal] /= a0[j+j*lda0];
                        }
                    }
                    else {
                        double scal = 1.0/a0[j+j*lda0];
                        if (l == 0)
                            cblas_dscal(mva0-j-1, (scal),
                                        &a0[j+1+j*lda0], 1);
                        else
                            cblas_dscal(mval, (scal),
                                        &al[j*ldal], 1);
                    }
                }

                // trailing update
                double zmone = -1.0;
                if (l == 0) {
                    cblas_dger(CblasColMajor,
                                mva0-j-1, k+kb-j-1,
                                (zmone), &a0[j+1+j*lda0], 1,
                                                    &a0[j+(j+1)*lda0], lda0,
                                                    &a0[j+1+(j+1)*lda0], lda0);
                }
                else {
                    cblas_dger(CblasColMajor,
                                mval, k+kb-j-1,
                                (zmone), &al[+j*ldal], 1,
                                                    &a0[j+(j+1)*lda0], lda0,
                                                    &al[+(j+1)*ldal], ldal);
                }
            }
            plasma_barrier_wait(barrier);
        }

        //===================================
        // right pivoting and trsm (rank 0)
        //===================================
        plasma_barrier_wait(barrier);
        if (rank == 0) {
            // pivot adjustment
            for (int i = k+1; i <= imin(A.m, k+kb); i++)
                ipiv[i-1] += k;

            // right pivoting
            for (int i = k; i < k+kb; i++) {
                double *ap = A((ipiv[i]-1)/A.mb, 0);
                int ldap = plasma_tile_mmain(A, (ipiv[i]-1)/A.mb);

                cblas_dswap(nva0-k-kb,
                            &a0[i+(k+kb)*lda0], lda0,
                            &ap[(ipiv[i]-1)%A.mb+(k+kb)*ldap], ldap);
            }
            // trsm
            double zone = 1.0;
            cblas_dtrsm(CblasColMajor,
                        PlasmaLeft, PlasmaLower,
                        PlasmaNoTrans, PlasmaUnit,
                        kb,
                        nva0-k-kb,
                        (zone), &a0[k+k*lda0], lda0,
                                           &a0[k+(k+kb)*lda0], lda0);
        }
        plasma_barrier_wait(barrier);

        //===================
        // gemm (all ranks)
        //===================
        double zone = 1.0;
        double zmone = -1.0;
        for (int i = rank; i < A.mt; i += size) {
            double *ai = A(i, 0);
            int mvai = plasma_tile_mview(A, i);
            int ldai = plasma_tile_mmain(A, i);

            if (i == 0) {
                cblas_dgemm(CblasColMajor,
                            CblasNoTrans, CblasNoTrans,
                            mva0-k-kb,
                            nva0-k-kb,
                            kb,
                            (zmone), &a0[k+kb+k*lda0], lda0,
                                                &a0[k+(k+kb)*lda0], lda0,
                            (zone),  &a0[(k+kb)+(k+kb)*lda0], lda0);
            }
            else {
                cblas_dgemm(CblasColMajor,
                            CblasNoTrans, CblasNoTrans,
                            mvai,
                            nva0-k-kb,
                            kb,
                            (zmone), &ai[k*ldai], ldai,
                                                &a0[k+(k+kb)*lda0], lda0,
                            (zone),  &ai[(k+kb)*ldai], ldai);
            }
        }
        plasma_barrier_wait(barrier);
    }

    //============================
    // left pivoting (all ranks)
    //============================
    for (int k = ib; k < imin(A.m, A.n); k += ib) {
        if (k%ib == rank) {
            for (int i = k; i < imin(A.m, A.n); i++) {
                double *ai = A(i/A.mb, 0);
                double *ap = A((ipiv[i]-1)/A.mb, 0);
                int ldai = plasma_tile_mmain(A, (i/A.mb));
                int ldap = plasma_tile_mmain(A, (ipiv[i]-1)/A.mb);

                cblas_dswap(ib,
                            &ai[i%A.mb+(k-ib)*ldai], ldai,
                            &ap[(ipiv[i]-1)%A.mb+(k-ib)*ldap], ldap);
            }
        }
    }

    // Free arrays for parallel max search.
    if (rank == 0) {
        free(max_idx);
        free(max_val);
    }

    // Only rank 0 returns errors.
    if (rank == 0)
        return info;
    else
        return 0;
}
