/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from core_blas/core_zgessq.c, normal z -> d, Thu Mar 10 18:58:58 2022
 *
 **/

#include <math.h>

#include "core_blas.h"
#include "plasma_types.h"
#include "core_lapack.h"

/******************************************************************************/
void core_dgessq(int m, int n,
                 const double *A, int lda,
                 double *scale, double *sumsq)
{
    int ione = 1;
    for (int j = 0; j < n; j++)
        // TODO: Inline this operation.
        LAPACK_dlassq(&m, &A[j*lda], &ione, scale, sumsq);
}

/******************************************************************************/
void core_omp_dgessq(int m, int n,
                     const double *A, int lda,
                     double *scale, double *sumsq,
                     plasma_sequence_t *sequence, plasma_request_t *request)
{
    #pragma omp task depend(in:A[0:lda*n]) \
                     depend(out:scale[0:n]) \
                     depend(out:sumsq[0:n])
    {
        if (sequence->status == PlasmaSuccess) {
            *scale = 0.0;
            *sumsq = 1.0;
            core_dgessq(m, n, A, lda, scale, sumsq);
        }
    }
}

/******************************************************************************/
void core_omp_dgessq_aux(int n,
                         const double *scale, const double *sumsq,
                         double *value,
                         plasma_sequence_t *sequence, plasma_request_t *request)
{
    #pragma omp task depend(in:scale[0:n]) \
                     depend(in:sumsq[0:n]) \
                     depend(out:value[0:1])
    {
        if (sequence->status == PlasmaSuccess) {
            double scl = 0.0;
            double sum = 1.0;
            for (int i = 0; i < n; i++) {
                if (scl < scale[i]) {
                    sum = sumsq[i] + sum*((scl/scale[i])*(scl/scale[i]));
                    scl = scale[i];
                }
                else {
                    sum = sum + sumsq[i]*(scale[i]/scl)*(scale[i]/scl);
                }
            }
            *value = scl*sqrt(sum);
        }
    }
}
