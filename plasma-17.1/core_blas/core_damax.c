/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from core_blas/core_dzamax.c, normal z -> d, Tue Feb  8 19:16:13 2022
 *
 **/

#include "core_blas.h"
#include "plasma_types.h"
#include "core_lapack.h"

#include <math.h>

/******************************************************************************/
void core_omp_damax(int colrow, int m, int n,
                     const double *A, int lda,
                     double *values,
                     plasma_sequence_t *sequence, plasma_request_t *request)
{
    switch (colrow) {
    case PlasmaColumnwise:
        #pragma omp task depend(in:A[0:lda*n]) \
                         depend(out:values[0:n])
        {
            if (sequence->status == PlasmaSuccess) {
                for (int j = 0; j < n; j++) {
                    values[j] = fabs(A[lda*j]);
                    for (int i = 1; i < m; i++) {
                        double tmp = fabs(A[lda*j+i]);
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
                    values[i] = fabs(A[i]);

                for (int j = 1; j < n; j++) {
                    for (int i = 0; i < m; i++) {
                        double tmp = fabs(A[lda*j+i]);
                        if (tmp > values[i])
                            values[i] = tmp;
                    }
                }
            }
        }
        break;
    }
}
