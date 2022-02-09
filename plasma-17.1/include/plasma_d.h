/**
 *
 * @file
 *
 *  PLASMA header.
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of Manchester, Univ. of California Berkeley and
 *  Univ. of Colorado Denver.
 *
 * @generated from include/plasma_z.h, normal z -> d, Tue Feb  8 19:15:03 2022
 *
 **/
#ifndef ICL_PLASMA_D_H
#define ICL_PLASMA_D_H

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
int plasma_damax(plasma_enum_t colrow,
                  int m, int n,
                  double *pA, int lda, double *values);

int plasma_dgbsv(int n, int kl, int ku, int nrhs,
                 double *pAB, int ldab, int *ipiv,
                 double *pB,  int ldb);

int plasma_dgbtrf(int m, int n, int kl, int ku,
                  double *pA, int lda, int *ipiv);

int plasma_dgbtrs(plasma_enum_t transa, int n, int kl, int ku, int nrhs,
                  double *pAB, int ldab,
                  int *ipiv,
                  double *pB,  int ldb);

int plasma_dgeadd(plasma_enum_t transa,
                  int m, int n,
                  double alpha, double *pA, int lda,
                  double beta,  double *pB, int ldb);

int plasma_dgelqf(int m, int n,
                  double *pA, int lda,
                  plasma_desc_t *T);

int plasma_dgelqs(int m, int n, int nrhs,
                  double *pA, int lda,
                  plasma_desc_t T,
                  double *pB, int ldb);

int plasma_dgels(plasma_enum_t trans,
                 int m, int n, int nrhs,
                 double *pA, int lda,
                 plasma_desc_t *T,
                 double *pB, int ldb);

int plasma_dgemm(plasma_enum_t transa, plasma_enum_t transb,
                 int m, int n, int k,
                 double alpha, double *pA, int lda,
                                           double *pB, int ldb,
                 double beta,  double *pC, int ldc);

int plasma_dgeqrf(int m, int n,
                  double *pA, int lda,
                  plasma_desc_t *T);

int plasma_dgeqrs(int m, int n, int nrhs,
                  double *pA, int lda,
                  plasma_desc_t T,
                  double *pB, int ldb);

int plasma_dgesv(int n, int nrhs,
                 double *pA, int lda, int *ipiv,
                 double *pB, int ldb);

int plasma_dgetrf(int m, int n,
                  double *pA, int lda, int *ipiv);

int plasma_dgetri(int n, double *pA, int lda, int *ipiv);

int plasma_dgetri_aux(int n, double *pA, int lda);

int plasma_dgetrs(int n, int nrhs,
                  double *pA, int lda, int *ipiv,
                  double *pB, int ldb);

int plasma_dsymm(plasma_enum_t side, plasma_enum_t uplo,
                 int m, int n,
                 double alpha, double *pA, int lda,
                                           double *pB, int ldb,
                 double beta,  double *pC, int ldc);

int plasma_dsyr2k(plasma_enum_t uplo, plasma_enum_t trans,
                  int n, int k,
                  double alpha, double *pA, int lda,
                                            double *pB, int ldb,
                  double beta,              double *pC, int ldc);

int plasma_dsyrk(plasma_enum_t uplo, plasma_enum_t trans,
                 int n, int k,
                 double alpha, double *pA, int lda,
                 double beta,  double *pC, int ldc);

int plasma_dlacpy(plasma_enum_t uplo,
                  int m, int n,
                  double *pA, int lda,
                  double *pB, int ldb);

double plasma_dlange(plasma_enum_t norm,
                     int m, int n,
                     double *pA, int lda);

double plasma_dlansy(plasma_enum_t norm, plasma_enum_t uplo,
                     int n,
                     double *pA, int lda);

double plasma_dlansy(plasma_enum_t norm, plasma_enum_t uplo,
                     int n,
                     double *pA, int lda);

double plasma_dlantr(plasma_enum_t norm, plasma_enum_t uplo, plasma_enum_t diag,
                     int m, int n,
                     double *pA, int lda);

int plasma_dlascl(plasma_enum_t uplo,
                  double cfrom, double cto,
                  int m, int n,
                  double *pA, int lda);

int plasma_dlaset(plasma_enum_t uplo,
                  int m, int n,
                  double alpha, double beta,
                  double *pA, int lda);

int plasma_dlaswp(plasma_enum_t colrow,
                  int m, int n,
                  double *pA, int lda,
                  int *ipiv, int incx);

int plasma_dlauum(plasma_enum_t uplo, int n,
                  double *pA, int lda);

int plasma_dpbsv(plasma_enum_t uplo,
                 int n, int kd, int nrhs,
                 double *pAB, int ldab,
                 double *pB,  int ldb);

int plasma_dpbtrf(plasma_enum_t uplo,
                  int n, int kd,
                  double *pAB, int ldab);

int plasma_dpbtrs(plasma_enum_t uplo,
                  int n, int kd, int nrhs,
                  double *pAB, int ldab,
                  double *pB,  int ldb);

int plasma_dposv(plasma_enum_t uplo,
                 int n, int nrhs,
                 double *pA, int lda,
                 double *pB, int ldb);

int plasma_dpotrf(plasma_enum_t uplo,
                  int n,
                  double *pA, int lda);

int plasma_dpotri(plasma_enum_t uplo,
                  int n,
                  double *pA, int lda);

int plasma_dpotrs(plasma_enum_t uplo,
                  int n, int nrhs,
                  double *pA, int lda,
                  double *pB, int ldb);

int plasma_dsymm(plasma_enum_t side, plasma_enum_t uplo,
                 int m, int n,
                 double alpha, double *pA, int lda,
                                           double *pB, int ldb,
                 double beta,  double *pC, int ldc);

int plasma_dsyr2k(plasma_enum_t uplo, plasma_enum_t trans,
                  int n, int k,
                  double alpha, double *pA, int lda,
                                            double *pB, int ldb,
                  double beta,  double *pC, int ldc);

int plasma_dsyrk(plasma_enum_t uplo, plasma_enum_t trans,
                 int n, int k,
                 double alpha, double *pA, int lda,
                 double beta,  double *pC, int ldc);

int plasma_dtradd(plasma_enum_t uplo, plasma_enum_t transa,
                  int m, int n,
                  double alpha, double *pA, int lda,
                  double beta,  double *pB, int ldb);

int plasma_dtrmm(plasma_enum_t side, plasma_enum_t uplo,
                 plasma_enum_t transa, plasma_enum_t diag,
                 int m, int n,
                 double alpha, double *pA, int lda,
                                           double *pB, int ldb);

int plasma_dtrsm(plasma_enum_t side, plasma_enum_t uplo,
                 plasma_enum_t transa, plasma_enum_t diag,
                 int m, int n,
                 double alpha, double *pA, int lda,
                                           double *pB, int ldb);

int plasma_dtrtri(plasma_enum_t uplo, plasma_enum_t diag,
                  int n, double *pA, int lda);

int plasma_dorglq(int m, int n, int k,
                  double *pA, int lda,
                  plasma_desc_t T,
                  double *pQ, int ldq);

int plasma_dorgqr(int m, int n, int k,
                  double *pA, int lda,
                  plasma_desc_t T,
                  double *pQ, int ldq);

int plasma_dormlq(plasma_enum_t side, plasma_enum_t trans,
                  int m, int n, int k,
                  double *pA, int lda,
                  plasma_desc_t T,
                  double *pC, int ldc);

int plasma_dormqr(plasma_enum_t side, plasma_enum_t trans,
                  int m, int n, int k,
                  double *pA, int lda,
                  plasma_desc_t T,
                  double *pC, int ldc);

/***************************************************************************//**
 *  Tile asynchronous interface.
 **/
void plasma_omp_damax(plasma_enum_t colrow, plasma_desc_t A,
                       double *work, double *values,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_dgbsv(plasma_desc_t AB, int *ipiv, plasma_desc_t B,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_dgbtrf(plasma_desc_t A, int *ipiv,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_dgbtrs(plasma_enum_t transa, plasma_desc_t AB, int *ipiv,
                       plasma_desc_t B,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_ddesc2ge(plasma_desc_t A,
                         double *pA, int lda,
                         plasma_sequence_t *sequence,
                         plasma_request_t *request);

void plasma_omp_ddesc2pb(plasma_desc_t A,
                         double *pA, int lda,
                         plasma_sequence_t *sequence,
                         plasma_request_t *request);

void plasma_omp_dge2desc(double *pA, int lda,
                         plasma_desc_t A,
                         plasma_sequence_t *sequence,
                         plasma_request_t *request);

void plasma_omp_dgeadd(plasma_enum_t transa,
                       double alpha, plasma_desc_t A,
                       double beta,  plasma_desc_t B,
                       plasma_sequence_t *sequence, plasma_request_t  *request);

void plasma_omp_dgelqf(plasma_desc_t A, plasma_desc_t T,
                       plasma_workspace_t work,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_dgelqs(plasma_desc_t A, plasma_desc_t T,
                       plasma_desc_t B, plasma_workspace_t work,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_dgels(plasma_enum_t trans,
                      plasma_desc_t A, plasma_desc_t T,
                      plasma_desc_t B, plasma_workspace_t work,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_dgemm(plasma_enum_t transa, plasma_enum_t transb,
                      double alpha, plasma_desc_t A,
                                                plasma_desc_t B,
                      double beta,  plasma_desc_t C,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_dgeqrf(plasma_desc_t A, plasma_desc_t T,
                       plasma_workspace_t work,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_dgeqrs(plasma_desc_t A, plasma_desc_t T,
                       plasma_desc_t B, plasma_workspace_t work,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_dgesv(plasma_desc_t A, int *ipiv,
                      plasma_desc_t B,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_dgetrf(plasma_desc_t A, int *ipiv,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_dgetri(plasma_desc_t A, int *ipiv, plasma_desc_t W,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_dgetri_aux(plasma_desc_t A, plasma_desc_t W,
                           plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_dgetrs(plasma_desc_t A, int *ipiv,
                       plasma_desc_t B,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_dsymm(plasma_enum_t side, plasma_enum_t uplo,
                      double alpha, plasma_desc_t A,
                                                plasma_desc_t B,
                      double beta,  plasma_desc_t C,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_dsyr2k(plasma_enum_t uplo, plasma_enum_t trans,
                       double alpha, plasma_desc_t A,
                                                 plasma_desc_t B,
                       double beta,              plasma_desc_t C,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_dsyrk(plasma_enum_t uplo, plasma_enum_t trans,
                      double alpha, plasma_desc_t A,
                      double beta,  plasma_desc_t C,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_dlacpy(plasma_enum_t uplo, plasma_desc_t A, plasma_desc_t B,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_dlange(plasma_enum_t norm, plasma_desc_t A,
                       double *work, double *value,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_dlansy(plasma_enum_t norm, plasma_enum_t uplo, plasma_desc_t A,
                       double *work, double *value,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_dlansy(plasma_enum_t norm, plasma_enum_t uplo, plasma_desc_t A,
                       double *work, double *value,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_dlantr(plasma_enum_t norm, plasma_enum_t uplo,
                       plasma_enum_t diag, plasma_desc_t A,
                       double *work, double *value,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_dlascl(plasma_enum_t uplo,
                       double cfrom, double cto,
                       plasma_desc_t A,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_dlaset(plasma_enum_t uplo,
                       double alpha, double beta,
                       plasma_desc_t A,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_dlaswp(plasma_enum_t colrow,
                       plasma_desc_t A,
                       int *ipiv, int incx,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_dlauum(plasma_enum_t uplo,
                       plasma_desc_t A,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_dpb2desc(double *pA, int lda,
                         plasma_desc_t A,
                         plasma_sequence_t *sequence,
                         plasma_request_t *request);

void plasma_omp_dpbsv(plasma_enum_t uplo, plasma_desc_t AB, plasma_desc_t B,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_dpbtrf(plasma_enum_t uplo, plasma_desc_t AB,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_dpbtrs(plasma_enum_t uplo, plasma_desc_t AB, plasma_desc_t B,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_dposv(plasma_enum_t uplo, plasma_desc_t A, plasma_desc_t B,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_dpotrf(plasma_enum_t uplo, plasma_desc_t A,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_dpotri(plasma_enum_t uplo, plasma_desc_t A,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_dpotrs(plasma_enum_t uplo, plasma_desc_t A, plasma_desc_t B,
                        plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_dsymm(plasma_enum_t side, plasma_enum_t uplo,
                      double alpha, plasma_desc_t A,
                                                plasma_desc_t B,
                      double beta,  plasma_desc_t C,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_dsyr2k(plasma_enum_t uplo, plasma_enum_t trans,
                       double alpha, plasma_desc_t A,
                                                 plasma_desc_t B,
                       double beta,  plasma_desc_t C,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_dsyrk(plasma_enum_t uplo, plasma_enum_t trans,
                      double alpha, plasma_desc_t A,
                      double beta,  plasma_desc_t C,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_dtradd(plasma_enum_t uplo, plasma_enum_t transa,
                       double alpha, plasma_desc_t A,
                       double beta,  plasma_desc_t B,
                       plasma_sequence_t *sequence, plasma_request_t  *request);

void plasma_omp_dtrmm(plasma_enum_t side, plasma_enum_t uplo,
                      plasma_enum_t transa, plasma_enum_t diag,
                      double alpha, plasma_desc_t A,
                                                plasma_desc_t B,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_dtrsm(plasma_enum_t side, plasma_enum_t uplo,
                      plasma_enum_t transa, plasma_enum_t diag,
                      double alpha, plasma_desc_t A,
                                                plasma_desc_t B,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_dtrtri(plasma_enum_t uplo, plasma_enum_t diag,
                       plasma_desc_t A,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_dorglq(plasma_desc_t A, plasma_desc_t T,
                       plasma_desc_t Q, plasma_workspace_t work,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_dorgqr(plasma_desc_t A, plasma_desc_t T,
                       plasma_desc_t Q, plasma_workspace_t work,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_dormlq(plasma_enum_t side, plasma_enum_t trans,
                       plasma_desc_t A, plasma_desc_t T,
                       plasma_desc_t C, plasma_workspace_t work,
                       plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_omp_dormqr(plasma_enum_t side, plasma_enum_t trans,
                       plasma_desc_t A, plasma_desc_t T,
                       plasma_desc_t C, plasma_workspace_t work,
                       plasma_sequence_t *sequence, plasma_request_t *request);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // ICL_PLASMA_D_H
