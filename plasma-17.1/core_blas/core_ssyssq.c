/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from core_blas/core_zsyssq.c, normal z -> s, Mon Nov 22 19:22:15 2021
 *
 **/

#include "core_blas.h"
#include "plasma_types.h"
#include "core_lapack.h"

#include <math.h>

/******************************************************************************/
void core_ssyssq(plasma_enum_t uplo,
                 int n,
                 const float *A, int lda,
                 float *scale, float *sumsq)
{
    int ione = 1;
    if (uplo == PlasmaUpper) {
        for (int j = 1; j < n; j++)
            // TODO: Inline this operation.
            LAPACK_slassq(&j, &A[lda*j], &ione, scale, sumsq);
    }
    else { // PlasmaLower
        for (int j = 0; j < n-1; j++) {
            int len = n-j-1;
            // TODO: Inline this operation.
            LAPACK_slassq(&len, &A[lda*j+j+1], &ione, scale, sumsq);
        }
    }
    *sumsq *= 2.0;
    for (int i = 0; i < n; i++) {
        // diagonal is complex, don't ignore complex part
        float absa = fabsf(A[lda*i+i]);
        if (absa != 0.0) { // != propagates nan
            if (*scale < absa) {
                *sumsq = 1.0 + *sumsq*((*scale/absa)*(*scale/absa));
                *scale = absa;
            }
            else {
                *sumsq = *sumsq + ((absa/(*scale))*(absa/(*scale)));
            }
        }
    }
}

/******************************************************************************/
void core_omp_ssyssq(plasma_enum_t uplo,
                     int n,
                     const float *A, int lda,
                     float *scale, float *sumsq,
                     plasma_sequence_t *sequence, plasma_request_t *request)
{
    #pragma omp task depend(in:A[0:lda*n]) \
                     depend(out:scale[0:n]) \
                     depend(out:sumsq[0:n])
    {
        if (sequence->status == PlasmaSuccess) {
            *scale = 0.0;
            *sumsq = 1.0;
            core_ssyssq(uplo, n, A, lda, scale, sumsq);
        }
    }
}

/******************************************************************************/
void core_omp_ssyssq_aux(int m, int n,
                         const float *scale, const float *sumsq,
                         float *value,
                         plasma_sequence_t *sequence, plasma_request_t *request)
{
    #pragma omp task depend(in:scale[0:n]) \
                     depend(in:sumsq[0:n]) \
                     depend(out:value[0:1])
    {
        if (sequence->status == PlasmaSuccess) {
            float scl = 0.0;
            float sum = 1.0;
            for (int j = 0; j < n; j++) {
                for (int i = j+1; i < n; i++) {
                    int idx = m*j+i;
                    if (scl < scale[idx]) {
                        sum = sumsq[idx] +
                            sum*((scl/scale[idx])*(scl/scale[idx]));
                        scl = scale[idx];
                    }
                    else {
                        sum = sum +
                            sumsq[idx]*((scale[idx]/scl)*(scale[idx]/scl));
                    }
                }
            }
            sum = 2.0*sum;
            for (int j = 0; j < n; j++) {
                int idx = m*j+j;
                if (scl < scale[idx]) {
                    sum = sumsq[idx] + sum*((scl/scale[idx])*(scl/scale[idx]));
                    scl = scale[idx];
                }
                else {
                    sum = sum + sumsq[idx]*((scale[idx]/scl)*(scale[idx]/scl));
                }
            }
            *value = scl*sqrtf(sum);
        }
    }
}
