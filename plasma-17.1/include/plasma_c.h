/**
 *
 * @file
 *
 *  PLASMA header.
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of Manchester, Univ. of California Berkeley and
 *  Univ. of Colorado Denver.
 *
 * @generated from include/plasma_z.h, normal z -> c, Mon Nov 22 19:09:34 2021
 *
 **/
#ifndef ICL_PLASMA_C_H
#define ICL_PLASMA_C_H

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
int plasma_scamax(plasma_enum_t colrow,
                  int m, int n,
                  plasma_complex32_t *pA, int lda, float *values);

int plasma_cgbsv(int n, int kl, int ku, int nrhs,
                 plasma_complex32_t *pAB, int ldab, int *ipiv,
                 plasma_complex32_t *pB,  int ldb);

int plasma_cgbtrf(int m, int n, int kl, int ku,
                  plasma_complex32_t *pA, int lda, int *ipiv);

int plasma_cgbtrs(plasma_enum_t transa, int n, int kl, int ku, int nrhs,
                  plasma_complex32_t *pAB, int ldab,
                  int *ipiv,
                  plasma_complex32_t *pB,  int ldb);

int plasma_cgeadd(plasma_enum_t transa,
                  int m, int n,
                  plasma_complex32_t alpha, plasma_complex32_t *pA, int lda,
                  plasma_complex32_t beta,  plasma_complex32_t *pB, int ldb);

int plasma_cgelqf(int m, int n,
                  plasma_complex32_t *pA, int lda,
                  plasma_desc_t *T);

int plasma_cgelqs(int m, int n, int nrhs,
                  plasma_complex32_t *pA, int lda,
                  plasma_desc_t T,
                  plasma_complex32_t *pB, int ldb);

int plasma_cgels(plasma_enum_t trans,
                 int m, int n, int nrhs,
                 plasma_complex32_t *pA, int lda,
                 plasma_desc_t *T,
                 plasma_complex32_t *pB, int ldb);

int plasma_cgemm(plasma_enum_t transa, plasma_enum_t transb,
                 int m, int n, int k,
                 plasma_complex32_t alpha, plasma_complex32_t *pA, int lda,
                                           plasma_complex32_t *pB, int ldb,
                 plasma_complex32_t beta,  plasma_complex32_t *pC, int ldc);

int plasma_cgeqrf(int m, int n,
                  plasma_complex32_t *pA, int lda,
                  plasma_desc_t *T);

int plasma_cgeqrs(int m, int n, int nrhs,
                  plasma_complex32_t *pA, int lda,
                  plasma_desc_t T,
                  plasma_complex32_t *pB, int ldb);

int plasma_cgesv(int n, int nrhs,
                 plasma_complex32_t *pA, int lda, int *ipiv,
                 plasma_complex32_t *pB, int ldb);

int plasma_cgetrf(int m, int n,
                  plasma_complex32_t *pA, int lda, int *ipiv);

int plasma_cgetri(int n, plasma_complex32_t *pA, int lda, int *ipiv);

int plasma_cgetri_aux(int n, plasma_complex32_t *pA, int lda);

int plasma_cgetrs(int n, int nrhs,
                  plasma_complex32_t *pA, int lda, int *ipiv,
                  plasma_complex32_t *pB, int ldb);

int plasma_chemm(plasma_enum_t side, plasma_enum_t uplo,
                 int m, int n,
                 plasma_complex32_t alpha, plasma_complex32_t *pA, int lda,
                                           plasma_complex32_t *pB, int ldb,
                 plasma_complex32_t beta,  plasma_complex32_t *pC, int ldc);

int plasma_cher2k(plasma_enum_t uplo, plasma_enum_t trans,
                  int n, int k,
                  plasma_complex32_t alpha, plasma_complex32_t *pA, int lda,
                                            plasma_complex32_t *pB, int ldb,
                  float beta,              plasma_complex32_t *pC, int ldc);

int plasma_cherk(plasma_enum_t uplo, plasma_enum_t trans,
                 int n, int k,
                 float alpha, plasma_complex32_t *pA, int lda,
                 float beta,  plasma_complex32_t *pC, int ldc);

int plasma_clacpy(plasma_enum_t uplo,
                  int m, int n,
                  plasma_complex32_t *pA, int lda,
                  plasma_complex32_t *pB, int ldb);

float plasma_clange(plasma_enum_t norm,
                     int m, int n,
                     plasma_complex32_t *pA, int lda);

float plasma_clanhe(plasma_enum_t norm, plasma_enum_t uplo,
                     int n,
                     plasma_complex32_t *pA, int lda);

float plasma_clansy(plasma_enum_t norm, plasma_enum_t uplo,
                     int n,
                     plasma_complex32_t *pA, int lda);

float plasma_clantr(plasma_enum_t norm, plasma_enum_t uplo, plasma_enum_t diag,
                     int m, int n,
                     plasma_complex32_t *pA, int lda);

int plasma_clascl(plasma_enum_t uplo,
                  float cfrom, float cto,
                  int m, int n,
                  plasma_complex32_t *pA, int lda);

int plasma_claset(plasma_enum_t uplo,
                  int m, int n,
                  plasma_complex32_t alpha, plasma_complex32_t beta,
                  plasma_complex32_t *pA, int lda);

int plasma_claswp(plasma_enum_t colrow,
                  int m, int n,
                  plasma_complex32_t *pA, int lda,
                  int *ipiv, int incx);

int plasma_clauum(plasma_enum_t uplo, int n,
                  plasma_complex32_t *pA, int lda);

int plasma_cpbsv(plasma_enum_t uplo,
                 int n, int kd, int nrhs,
                 plasma_complex32_t *pAB, int ldab,
                 plasma_complex32_t *pB,  int ldb);

int plasma_cpbtrf(plasma_enum_t uplo,
                  int n, int kd,
                  plasma_complex32_t *pAB, int ldab);

int plasma_cpbtrs(plasma_enum_t uplo,
                  int n, int kd, int nrhs,
                  plasma_complex32_t *pAB, int ldab,
                  plasma_complex32_t *pB,  int ldb);

int plasma_cposv(plasma_enum_t uplo,
                 int n, int nrhs,
                 plasma_complex32_t *pA, int lda,
                 plasma_complex32_t *pB, int ldb);

int plasma_cpotrf(plasma_enum_t uplo,
                  int n,
                  plasma_complex32_t *pA, int lda);

int plasma_cpotri(plasma_enum_t uplo,
                  int n,
                  plasma_complex32_t *pA, int lda);

int plasma_cpotrs(plasma_enum_t uplo,
                  int n, int nrhs,
                  plasma_complex32_t *pA, int lda,
                  plasma_complex32_t *pB, int ldb);

int plasma_csymm(plasma_enum_t side, plasma_enum_t uplo,
                 int m, int n,
                 plasma_complex32_t alpha, plasma_complex32_t *pA, int lda,
                                           plasma_complex32_t *pB, int ldb,
                 plasma_complex32_t beta,  plasma_complex32_t *pC, int ldc);

int plasma_csyr2k(plasma_enum_t uplo, plasma_enum_t trans,
                  int n, int k,
                  plasma_complex32_t alpha, plasma_complex32_t *pA, int lda,
                                            plasma_complex32_t *pB, int ldb,
                  plasma_complex32_t beta,  plasma_complex32_t *pC, int ldc);

int plasma_csyrk(plasma_enum_t uplo, plasma_enum_t trans,
                 int n, int k,
                 plasma_complex32_t alpha, plasma_complex32_t *pA, int lda,
                 plasma_complex32_t beta,  plasma_complex32_t *pC, int ldc);

int plasma_ctradd(plasma_enum_t uplo, plasma_enum_t transa,
                  int m, int n,
                  plasma_complex32_t alpha, plasma_complex32_t *pA, int lda,
                  plasma_complex32_t beta,  plasma_complex32_t *pB, int ldb);

int plasma_ctrmm(plasma_enum_t side, plasma_enum_t uplo,
                 plasma_enum_t transa, plasma_enum_t diag,
                 int m, int n,
                 plasma_complex32_t alpha, plasma_complex32_t *pA, int lda,
                                           plasma_complex32_t *pB, int ldb);

int plasma_ctrsm(plasma_enum_t side, plasma_enum_t uplo,
                 plasma_enum_t transa, plasma_enum_t diag,
                 int m, int n,
                 plasma_complex32_t alpha, plasma_complex32_t *pA, int lda,
                                           plasma_complex32_t *pB, int ldb);

int plasma_ctrtri(plasma_enum_t uplo, plasma_enum_t diag,
                  int n, plasma_complex32_t *pA, int lda);

int plasma_cunglq(int m, int n, int k,
                  plasma_complex32_t *pA, int lda,
                  plasma_desc_t T,
                  plasma_complex32_t *pQ, int ldq);

int plasma_cungqr(int m, int n, int k,
                  plasma_complex32_t *pA, int lda,
                  plasma_desc_t T,
                  plasma_complex32_t *pQ, int ldq);

int plasma_cunmlq(plasma_enum_t side, plasma_enum_t trans,
                  int m, int n, int k,
                  plasma_complex32_t *pA, int lda,
                  plasma_desc_t T,
                  plasma_complex32_t *pC, int ldc);

int plasma_cunmqr(plasma_enum_t side, plasma_enum_t trans,
                  int m, int n, int k,
                  plasma_complex32_t *pA, int lda,
                  plasma_desc_t T,
                  plasma_complex32_t *pC, int ldc);

/***************************************************************************//**
 *  Tile asynchronous interface.
 **/
void plasma_omp_scamax(plasma_enum_t colrow, plasma_desc_t A,
                       float *work, float *values,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_cgbsv(plasma_desc_t AB, int *ipiv, plasma_desc_t B,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_cgbtrf(plasma_desc_t A, int *ipiv,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_cgbtrs(plasma_enum_t transa, plasma_desc_t AB, int *ipiv,
                       plasma_desc_t B,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_cdesc2ge(plasma_desc_t A,
                         plasma_complex32_t *pA, int lda,
                         plasma_sequence_t *sequence,
                         plasma_request_t *request);

void plasma_omp_cdesc2pb(plasma_desc_t A,
                         plasma_complex32_t *pA, int lda,
                         plasma_sequence_t *sequence,
                         plasma_request_t *request);

void plasma_omp_cge2desc(plasma_complex32_t *pA, int lda,
                         plasma_desc_t A,
                         plasma_sequence_t *sequence,
                         plasma_request_t *request);

void plasma_omp_cgeadd(plasma_enum_t transa,
                       plasma_complex32_t alpha, plasma_desc_t A,
                       plasma_complex32_t beta,  plasma_desc_t B,
                       plasma_sequence_t *sequence, plasma_request_t  *request);

void plasma_omp_cgelqf(plasma_desc_t A, plasma_desc_t T,
                       plasma_workspace_t work,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_cgelqs(plasma_desc_t A, plasma_desc_t T,
                       plasma_desc_t B, plasma_workspace_t work,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_cgels(plasma_enum_t trans,
                      plasma_desc_t A, plasma_desc_t T,
                      plasma_desc_t B, plasma_workspace_t work,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_cgemm(plasma_enum_t transa, plasma_enum_t transb,
                      plasma_complex32_t alpha, plasma_desc_t A,
                                                plasma_desc_t B,
                      plasma_complex32_t beta,  plasma_desc_t C,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_cgeqrf(plasma_desc_t A, plasma_desc_t T,
                       plasma_workspace_t work,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_cgeqrs(plasma_desc_t A, plasma_desc_t T,
                       plasma_desc_t B, plasma_workspace_t work,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_cgesv(plasma_desc_t A, int *ipiv,
                      plasma_desc_t B,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_cgetrf(plasma_desc_t A, int *ipiv,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_cgetri(plasma_desc_t A, int *ipiv, plasma_desc_t W,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_cgetri_aux(plasma_desc_t A, plasma_desc_t W,
                           plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_cgetrs(plasma_desc_t A, int *ipiv,
                       plasma_desc_t B,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_chemm(plasma_enum_t side, plasma_enum_t uplo,
                      plasma_complex32_t alpha, plasma_desc_t A,
                                                plasma_desc_t B,
                      plasma_complex32_t beta,  plasma_desc_t C,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_cher2k(plasma_enum_t uplo, plasma_enum_t trans,
                       plasma_complex32_t alpha, plasma_desc_t A,
                                                 plasma_desc_t B,
                       float beta,              plasma_desc_t C,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_cherk(plasma_enum_t uplo, plasma_enum_t trans,
                      float alpha, plasma_desc_t A,
                      float beta,  plasma_desc_t C,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_clacpy(plasma_enum_t uplo, plasma_desc_t A, plasma_desc_t B,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_clange(plasma_enum_t norm, plasma_desc_t A,
                       float *work, float *value,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_clanhe(plasma_enum_t norm, plasma_enum_t uplo, plasma_desc_t A,
                       float *work, float *value,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_clansy(plasma_enum_t norm, plasma_enum_t uplo, plasma_desc_t A,
                       float *work, float *value,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_clantr(plasma_enum_t norm, plasma_enum_t uplo,
                       plasma_enum_t diag, plasma_desc_t A,
                       float *work, float *value,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_clascl(plasma_enum_t uplo,
                       float cfrom, float cto,
                       plasma_desc_t A,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_claset(plasma_enum_t uplo,
                       plasma_complex32_t alpha, plasma_complex32_t beta,
                       plasma_desc_t A,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_claswp(plasma_enum_t colrow,
                       plasma_desc_t A,
                       int *ipiv, int incx,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_clauum(plasma_enum_t uplo,
                       plasma_desc_t A,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_cpb2desc(plasma_complex32_t *pA, int lda,
                         plasma_desc_t A,
                         plasma_sequence_t *sequence,
                         plasma_request_t *request);

void plasma_omp_cpbsv(plasma_enum_t uplo, plasma_desc_t AB, plasma_desc_t B,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_cpbtrf(plasma_enum_t uplo, plasma_desc_t AB,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_cpbtrs(plasma_enum_t uplo, plasma_desc_t AB, plasma_desc_t B,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_cposv(plasma_enum_t uplo, plasma_desc_t A, plasma_desc_t B,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_cpotrf(plasma_enum_t uplo, plasma_desc_t A,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_cpotri(plasma_enum_t uplo, plasma_desc_t A,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_cpotrs(plasma_enum_t uplo, plasma_desc_t A, plasma_desc_t B,
                        plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_csymm(plasma_enum_t side, plasma_enum_t uplo,
                      plasma_complex32_t alpha, plasma_desc_t A,
                                                plasma_desc_t B,
                      plasma_complex32_t beta,  plasma_desc_t C,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_csyr2k(plasma_enum_t uplo, plasma_enum_t trans,
                       plasma_complex32_t alpha, plasma_desc_t A,
                                                 plasma_desc_t B,
                       plasma_complex32_t beta,  plasma_desc_t C,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_csyrk(plasma_enum_t uplo, plasma_enum_t trans,
                      plasma_complex32_t alpha, plasma_desc_t A,
                      plasma_complex32_t beta,  plasma_desc_t C,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_ctradd(plasma_enum_t uplo, plasma_enum_t transa,
                       plasma_complex32_t alpha, plasma_desc_t A,
                       plasma_complex32_t beta,  plasma_desc_t B,
                       plasma_sequence_t *sequence, plasma_request_t  *request);

void plasma_omp_ctrmm(plasma_enum_t side, plasma_enum_t uplo,
                      plasma_enum_t transa, plasma_enum_t diag,
                      plasma_complex32_t alpha, plasma_desc_t A,
                                                plasma_desc_t B,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_ctrsm(plasma_enum_t side, plasma_enum_t uplo,
                      plasma_enum_t transa, plasma_enum_t diag,
                      plasma_complex32_t alpha, plasma_desc_t A,
                                                plasma_desc_t B,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_ctrtri(plasma_enum_t uplo, plasma_enum_t diag,
                       plasma_desc_t A,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_cunglq(plasma_desc_t A, plasma_desc_t T,
                       plasma_desc_t Q, plasma_workspace_t work,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_cungqr(plasma_desc_t A, plasma_desc_t T,
                       plasma_desc_t Q, plasma_workspace_t work,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_cunmlq(plasma_enum_t side, plasma_enum_t trans,
                       plasma_desc_t A, plasma_desc_t T,
                       plasma_desc_t C, plasma_workspace_t work,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_cunmqr(plasma_enum_t side, plasma_enum_t trans,
                       plasma_desc_t A, plasma_desc_t T,
                       plasma_desc_t C, plasma_workspace_t work,
                       plasma_sequence_t *sequence, plasma_request_t *request);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // ICL_PLASMA_C_H
