/**
 *
 * @file
 *
 *  PLASMA header.
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of Manchester, Univ. of California Berkeley and
 *  Univ. of Colorado Denver.
 *
 * @generated from include/plasma_z.h, normal z -> s, Mon Nov 22 19:09:34 2021
 *
 **/
#ifndef ICL_PLASMA_S_H
#define ICL_PLASMA_S_H

#include "plasma_async.h"
#include "plasma_barrier.h"
#include "plasma_descriptor.h"
#include "plasma_workspace.h"

#ifdef __cplusplus
extern "C" {
#endif

/***************************************************************************//**
 *  Standard interface.
 **/
int plasma_samax(plasma_enum_t colrow,
                  int m, int n,
                  float *pA, int lda, float *values);

int plasma_sgbsv(int n, int kl, int ku, int nrhs,
                 float *pAB, int ldab, int *ipiv,
                 float *pB,  int ldb);

int plasma_sgbtrf(int m, int n, int kl, int ku,
                  float *pA, int lda, int *ipiv);

int plasma_sgbtrs(plasma_enum_t transa, int n, int kl, int ku, int nrhs,
                  float *pAB, int ldab,
                  int *ipiv,
                  float *pB,  int ldb);

int plasma_sgeadd(plasma_enum_t transa,
                  int m, int n,
                  float alpha, float *pA, int lda,
                  float beta,  float *pB, int ldb);

int plasma_sgelqf(int m, int n,
                  float *pA, int lda,
                  plasma_desc_t *T);

int plasma_sgelqs(int m, int n, int nrhs,
                  float *pA, int lda,
                  plasma_desc_t T,
                  float *pB, int ldb);

int plasma_sgels(plasma_enum_t trans,
                 int m, int n, int nrhs,
                 float *pA, int lda,
                 plasma_desc_t *T,
                 float *pB, int ldb);

int plasma_sgemm(plasma_enum_t transa, plasma_enum_t transb,
                 int m, int n, int k,
                 float alpha, float *pA, int lda,
                                           float *pB, int ldb,
                 float beta,  float *pC, int ldc);

int plasma_sgeqrf(int m, int n,
                  float *pA, int lda,
                  plasma_desc_t *T);

int plasma_sgeqrs(int m, int n, int nrhs,
                  float *pA, int lda,
                  plasma_desc_t T,
                  float *pB, int ldb);

int plasma_sgesv(int n, int nrhs,
                 float *pA, int lda, int *ipiv,
                 float *pB, int ldb);

int plasma_sgetrf(int m, int n,
                  float *pA, int lda, int *ipiv);

int plasma_sgetri(int n, float *pA, int lda, int *ipiv);

int plasma_sgetri_aux(int n, float *pA, int lda);

int plasma_sgetrs(int n, int nrhs,
                  float *pA, int lda, int *ipiv,
                  float *pB, int ldb);

int plasma_ssymm(plasma_enum_t side, plasma_enum_t uplo,
                 int m, int n,
                 float alpha, float *pA, int lda,
                                           float *pB, int ldb,
                 float beta,  float *pC, int ldc);

int plasma_ssyr2k(plasma_enum_t uplo, plasma_enum_t trans,
                  int n, int k,
                  float alpha, float *pA, int lda,
                                            float *pB, int ldb,
                  float beta,              float *pC, int ldc);

int plasma_ssyrk(plasma_enum_t uplo, plasma_enum_t trans,
                 int n, int k,
                 float alpha, float *pA, int lda,
                 float beta,  float *pC, int ldc);

int plasma_slacpy(plasma_enum_t uplo,
                  int m, int n,
                  float *pA, int lda,
                  float *pB, int ldb);

float plasma_slange(plasma_enum_t norm,
                     int m, int n,
                     float *pA, int lda);

float plasma_slansy(plasma_enum_t norm, plasma_enum_t uplo,
                     int n,
                     float *pA, int lda);

float plasma_slansy(plasma_enum_t norm, plasma_enum_t uplo,
                     int n,
                     float *pA, int lda);

float plasma_slantr(plasma_enum_t norm, plasma_enum_t uplo, plasma_enum_t diag,
                     int m, int n,
                     float *pA, int lda);

int plasma_slascl(plasma_enum_t uplo,
                  float cfrom, float cto,
                  int m, int n,
                  float *pA, int lda);

int plasma_slaset(plasma_enum_t uplo,
                  int m, int n,
                  float alpha, float beta,
                  float *pA, int lda);

int plasma_slaswp(plasma_enum_t colrow,
                  int m, int n,
                  float *pA, int lda,
                  int *ipiv, int incx);

int plasma_slauum(plasma_enum_t uplo, int n,
                  float *pA, int lda);

int plasma_spbsv(plasma_enum_t uplo,
                 int n, int kd, int nrhs,
                 float *pAB, int ldab,
                 float *pB,  int ldb);

int plasma_spbtrf(plasma_enum_t uplo,
                  int n, int kd,
                  float *pAB, int ldab);

int plasma_spbtrs(plasma_enum_t uplo,
                  int n, int kd, int nrhs,
                  float *pAB, int ldab,
                  float *pB,  int ldb);

int plasma_sposv(plasma_enum_t uplo,
                 int n, int nrhs,
                 float *pA, int lda,
                 float *pB, int ldb);

int plasma_spotrf(plasma_enum_t uplo,
                  int n,
                  float *pA, int lda);

int plasma_spotri(plasma_enum_t uplo,
                  int n,
                  float *pA, int lda);

int plasma_spotrs(plasma_enum_t uplo,
                  int n, int nrhs,
                  float *pA, int lda,
                  float *pB, int ldb);

int plasma_ssymm(plasma_enum_t side, plasma_enum_t uplo,
                 int m, int n,
                 float alpha, float *pA, int lda,
                                           float *pB, int ldb,
                 float beta,  float *pC, int ldc);

int plasma_ssyr2k(plasma_enum_t uplo, plasma_enum_t trans,
                  int n, int k,
                  float alpha, float *pA, int lda,
                                            float *pB, int ldb,
                  float beta,  float *pC, int ldc);

int plasma_ssyrk(plasma_enum_t uplo, plasma_enum_t trans,
                 int n, int k,
                 float alpha, float *pA, int lda,
                 float beta,  float *pC, int ldc);

int plasma_stradd(plasma_enum_t uplo, plasma_enum_t transa,
                  int m, int n,
                  float alpha, float *pA, int lda,
                  float beta,  float *pB, int ldb);

int plasma_strmm(plasma_enum_t side, plasma_enum_t uplo,
                 plasma_enum_t transa, plasma_enum_t diag,
                 int m, int n,
                 float alpha, float *pA, int lda,
                                           float *pB, int ldb);

int plasma_strsm(plasma_enum_t side, plasma_enum_t uplo,
                 plasma_enum_t transa, plasma_enum_t diag,
                 int m, int n,
                 float alpha, float *pA, int lda,
                                           float *pB, int ldb);

int plasma_strtri(plasma_enum_t uplo, plasma_enum_t diag,
                  int n, float *pA, int lda);

int plasma_sorglq(int m, int n, int k,
                  float *pA, int lda,
                  plasma_desc_t T,
                  float *pQ, int ldq);

int plasma_sorgqr(int m, int n, int k,
                  float *pA, int lda,
                  plasma_desc_t T,
                  float *pQ, int ldq);

int plasma_sormlq(plasma_enum_t side, plasma_enum_t trans,
                  int m, int n, int k,
                  float *pA, int lda,
                  plasma_desc_t T,
                  float *pC, int ldc);

int plasma_sormqr(plasma_enum_t side, plasma_enum_t trans,
                  int m, int n, int k,
                  float *pA, int lda,
                  plasma_desc_t T,
                  float *pC, int ldc);

/***************************************************************************//**
 *  Tile asynchronous interface.
 **/
void plasma_omp_samax(plasma_enum_t colrow, plasma_desc_t A,
                       float *work, float *values,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_sgbsv(plasma_desc_t AB, int *ipiv, plasma_desc_t B,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_sgbtrf(plasma_desc_t A, int *ipiv,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_sgbtrs(plasma_enum_t transa, plasma_desc_t AB, int *ipiv,
                       plasma_desc_t B,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_sdesc2ge(plasma_desc_t A,
                         float *pA, int lda,
                         plasma_sequence_t *sequence,
                         plasma_request_t *request);

void plasma_omp_sdesc2pb(plasma_desc_t A,
                         float *pA, int lda,
                         plasma_sequence_t *sequence,
                         plasma_request_t *request);

void plasma_omp_sge2desc(float *pA, int lda,
                         plasma_desc_t A,
                         plasma_sequence_t *sequence,
                         plasma_request_t *request);

void plasma_omp_sgeadd(plasma_enum_t transa,
                       float alpha, plasma_desc_t A,
                       float beta,  plasma_desc_t B,
                       plasma_sequence_t *sequence, plasma_request_t  *request);

void plasma_omp_sgelqf(plasma_desc_t A, plasma_desc_t T,
                       plasma_workspace_t work,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_sgelqs(plasma_desc_t A, plasma_desc_t T,
                       plasma_desc_t B, plasma_workspace_t work,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_sgels(plasma_enum_t trans,
                      plasma_desc_t A, plasma_desc_t T,
                      plasma_desc_t B, plasma_workspace_t work,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_sgemm(plasma_enum_t transa, plasma_enum_t transb,
                      float alpha, plasma_desc_t A,
                                                plasma_desc_t B,
                      float beta,  plasma_desc_t C,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_sgeqrf(plasma_desc_t A, plasma_desc_t T,
                       plasma_workspace_t work,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_sgeqrs(plasma_desc_t A, plasma_desc_t T,
                       plasma_desc_t B, plasma_workspace_t work,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_sgesv(plasma_desc_t A, int *ipiv,
                      plasma_desc_t B,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_sgetrf(plasma_desc_t A, int *ipiv,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_sgetri(plasma_desc_t A, int *ipiv, plasma_desc_t W,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_sgetri_aux(plasma_desc_t A, plasma_desc_t W,
                           plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_sgetrs(plasma_desc_t A, int *ipiv,
                       plasma_desc_t B,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_ssymm(plasma_enum_t side, plasma_enum_t uplo,
                      float alpha, plasma_desc_t A,
                                                plasma_desc_t B,
                      float beta,  plasma_desc_t C,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_ssyr2k(plasma_enum_t uplo, plasma_enum_t trans,
                       float alpha, plasma_desc_t A,
                                                 plasma_desc_t B,
                       float beta,              plasma_desc_t C,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_ssyrk(plasma_enum_t uplo, plasma_enum_t trans,
                      float alpha, plasma_desc_t A,
                      float beta,  plasma_desc_t C,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_slacpy(plasma_enum_t uplo, plasma_desc_t A, plasma_desc_t B,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_slange(plasma_enum_t norm, plasma_desc_t A,
                       float *work, float *value,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_slansy(plasma_enum_t norm, plasma_enum_t uplo, plasma_desc_t A,
                       float *work, float *value,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_slansy(plasma_enum_t norm, plasma_enum_t uplo, plasma_desc_t A,
                       float *work, float *value,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_slantr(plasma_enum_t norm, plasma_enum_t uplo,
                       plasma_enum_t diag, plasma_desc_t A,
                       float *work, float *value,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_slascl(plasma_enum_t uplo,
                       float cfrom, float cto,
                       plasma_desc_t A,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_slaset(plasma_enum_t uplo,
                       float alpha, float beta,
                       plasma_desc_t A,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_slaswp(plasma_enum_t colrow,
                       plasma_desc_t A,
                       int *ipiv, int incx,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_slauum(plasma_enum_t uplo,
                       plasma_desc_t A,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_spb2desc(float *pA, int lda,
                         plasma_desc_t A,
                         plasma_sequence_t *sequence,
                         plasma_request_t *request);

void plasma_omp_spbsv(plasma_enum_t uplo, plasma_desc_t AB, plasma_desc_t B,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_spbtrf(plasma_enum_t uplo, plasma_desc_t AB,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_spbtrs(plasma_enum_t uplo, plasma_desc_t AB, plasma_desc_t B,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_sposv(plasma_enum_t uplo, plasma_desc_t A, plasma_desc_t B,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_spotrf(plasma_enum_t uplo, plasma_desc_t A,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_spotri(plasma_enum_t uplo, plasma_desc_t A,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_spotrs(plasma_enum_t uplo, plasma_desc_t A, plasma_desc_t B,
                        plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_ssymm(plasma_enum_t side, plasma_enum_t uplo,
                      float alpha, plasma_desc_t A,
                                                plasma_desc_t B,
                      float beta,  plasma_desc_t C,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_ssyr2k(plasma_enum_t uplo, plasma_enum_t trans,
                       float alpha, plasma_desc_t A,
                                                 plasma_desc_t B,
                       float beta,  plasma_desc_t C,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_ssyrk(plasma_enum_t uplo, plasma_enum_t trans,
                      float alpha, plasma_desc_t A,
                      float beta,  plasma_desc_t C,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_stradd(plasma_enum_t uplo, plasma_enum_t transa,
                       float alpha, plasma_desc_t A,
                       float beta,  plasma_desc_t B,
                       plasma_sequence_t *sequence, plasma_request_t  *request);

void plasma_omp_strmm(plasma_enum_t side, plasma_enum_t uplo,
                      plasma_enum_t transa, plasma_enum_t diag,
                      float alpha, plasma_desc_t A,
                                                plasma_desc_t B,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_strsm(plasma_enum_t side, plasma_enum_t uplo,
                      plasma_enum_t transa, plasma_enum_t diag,
                      float alpha, plasma_desc_t A,
                                                plasma_desc_t B,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_strtri(plasma_enum_t uplo, plasma_enum_t diag,
                       plasma_desc_t A,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_sorglq(plasma_desc_t A, plasma_desc_t T,
                       plasma_desc_t Q, plasma_workspace_t work,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_sorgqr(plasma_desc_t A, plasma_desc_t T,
                       plasma_desc_t Q, plasma_workspace_t work,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_sormlq(plasma_enum_t side, plasma_enum_t trans,
                       plasma_desc_t A, plasma_desc_t T,
                       plasma_desc_t C, plasma_workspace_t work,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_sormqr(plasma_enum_t side, plasma_enum_t trans,
                       plasma_desc_t A, plasma_desc_t T,
                       plasma_desc_t C, plasma_workspace_t work,
                       plasma_sequence_t *sequence, plasma_request_t *request);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // ICL_PLASMA_S_H
