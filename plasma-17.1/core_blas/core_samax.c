/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from core_blas/core_dzamax.c, normal z -> s, Mon Nov 22 19:22:27 2021
 *
 **/

#include "core_blas.h"
#include "plasma_types.h"
#include "core_lapack.h"

#include <math.h>

/******************************************************************************/
void core_omp_samax(int colrow, int m, int n,
                     const float *A, int lda,
                     float *values,
                     plasma_sequence_t *sequence, plasma_request_t *request)
{
    switch (colrow) {
    case PlasmaColumnwise:
        #pragma omp task depend(in:A[0:lda*n]) \
                         depend(out:values[0:n])
        {
            if (sequence->status == PlasmaSuccess) {
                for (int j = 0; j < n; j++) {
                    values[j] = fabsf(A[lda*j]);
                    for (int i = 1; i < m; i++) {
                        float tmp = fabsf(A[lda*j+i]);
                        if (tmp > values[j])
                            values[j] = tmp;
                    }
                }
            }
        }
        break;
    case PlasmaRowwise:
        #pragma omp task depend(in:A[0:lda*n]) \
                         depend(out:values[0:m])
        {
            if (sequence->status == PlasmaSuccess) {
                for (int i = 0; i < m; i++)
                    values[i] = fabsf(A[i]);

                for (int j = 1; j < n; j++) {
                    for (int i = 0; i < m; i++) {
                        float tmp = fabsf(A[lda*j+i]);
                        if (tmp > values[i])
                            values[i] = tmp;
                    }
                }
            }
        }
        break;
    }
}
