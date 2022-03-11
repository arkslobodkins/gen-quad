/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from include/core_blas_zc.h, mixed zc -> ds, Thu Mar 10 16:29:49 2022
 *
 **/
#ifndef ICL_CORE_BLAS_DS_H
#define ICL_CORE_BLAS_DS_H

#include "plasma_async.h"
#include "plasma_types.h"
#include "plasma_workspace.h"

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************/
void core_dlag2s(int m, int n,
                 double *A,  int lda,
                 float *As, int ldas);

void core_slag2d(int m, int n,
                 float *As, int ldas,
                 double *A,  int lda);

/******************************************************************************/
void core_omp_dlag2s(int m, int n,
                     double *A,  int lda,
                     float *As, int ldas,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_slag2d(int m, int n,
                     float *As, int ldas,
                     double *A,  int lda,
                     plasma_sequence_t *sequence, plasma_request_t *request);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // ICL_CORE_BLAS_DS_H
