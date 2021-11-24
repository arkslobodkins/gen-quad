/**
 *
 * @file
 *
 *  PLASMA header.
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of Manchester, Univ. of California Berkeley and
 *  Univ. of Colorado Denver.
 *
 * @generated from include/plasma_zc.h, mixed zc -> ds, Mon Nov 22 19:09:34 2021
 *
 **/
#ifndef ICL_PLASMA_DS_H
#define ICL_PLASMA_DS_H

#include "plasma_async.h"
#include "plasma_descriptor.h"
#include "plasma_workspace.h"

#ifdef __cplusplus
extern "C" {
#endif

/***************************************************************************//**
 *  Standard interface
 **/
int plasma_dsgesv(int n, int nrhs,
                  double *pA, int lda, int *ipiv,
                  double *pB, int ldb,
                  double *pX, int ldx, int *iter);

int plasma_dsposv(plasma_enum_t uplo, int n, int nrhs,
                  double *pA, int lda,
                  double *pB, int ldb,
                  double *pX, int ldx, int *iter);

int plasma_dlag2s(int m, int n,
                  double *pA,  int lda,
                  float *pAs, int ldas);

int plasma_slag2d(int m, int n,
                  float *pAs, int ldas,
                  double *pA,  int lda);

/***************************************************************************//**
 *  Tile asynchronous interface
 **/
void plasma_omp_dsgesv(plasma_desc_t A,  int *ipiv,
                       plasma_desc_t B,  plasma_desc_t X,
                       plasma_desc_t As, plasma_desc_t Xs, plasma_desc_t R,
                       double *work, double *Rnorm, double *Xnorm, int *iter,
                       plasma_sequence_t *sequence,
                       plasma_request_t  *request);

void plasma_omp_dsposv(plasma_enum_t uplo,
                       plasma_desc_t A,  plasma_desc_t B,  plasma_desc_t X,
                       plasma_desc_t As, plasma_desc_t Xs, plasma_desc_t R,
                       double *W,  double *Rnorm, double *Xnorm, int *iter,
                       plasma_sequence_t *sequence,
                       plasma_request_t  *request);

void plasma_omp_dlag2s(plasma_desc_t A, plasma_desc_t As,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_slag2d(plasma_desc_t As, plasma_desc_t A,
                       plasma_sequence_t *sequence, plasma_request_t *request);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // ICL_PLASMA_DS_H
