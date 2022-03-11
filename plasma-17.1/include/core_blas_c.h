/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from include/core_blas_z.h, normal z -> c, Thu Mar 10 16:29:49 2022
 *
 **/
#ifndef ICL_CORE_BLAS_C_H
#define ICL_CORE_BLAS_C_H

#include "plasma_async.h"
#include "plasma_barrier.h"
#include "plasma_descriptor.h"
#include "plasma_types.h"
#include "plasma_workspace.h"
#include "plasma_descriptor.h"

#ifdef __cplusplus
extern "C" {
#endif

#define COMPLEX

/******************************************************************************/
#ifdef COMPLEX
float core_scabs1(plasma_complex32_t alpha);
#endif

int core_cgeadd(plasma_enum_t transa,
                int m, int n,
                plasma_complex32_t alpha, const plasma_complex32_t *A, int lda,
                plasma_complex32_t beta,        plasma_complex32_t *B, int ldb);

int core_cgelqt(int m, int n, int ib,
                plasma_complex32_t *A, int lda,
                plasma_complex32_t *T, int ldt,
                plasma_complex32_t *tau,
                plasma_complex32_t *work);

void core_cgemm(plasma_enum_t transa, plasma_enum_t transb,
                int m, int n, int k,
                plasma_complex32_t alpha, const plasma_complex32_t *A, int lda,
                                          const plasma_complex32_t *B, int ldb,
                plasma_complex32_t beta,        plasma_complex32_t *C, int ldc);

int core_cgeqrt(int m, int n, int ib,
                plasma_complex32_t *A, int lda,
                plasma_complex32_t *T, int ldt,
                plasma_complex32_t *tau,
                plasma_complex32_t *work);

void core_cgessq(int m, int n,
                 const plasma_complex32_t *A, int lda,
                 float *scale, float *sumsq);

int core_cgetrf(plasma_desc_t A, int *ipiv, int ib, int rank, int size,
                plasma_barrier_t *barrier);

void core_chemm(plasma_enum_t side, plasma_enum_t uplo,
                int m, int n,
                plasma_complex32_t alpha, const plasma_complex32_t *A, int lda,
                                          const plasma_complex32_t *B, int ldb,
                plasma_complex32_t beta,        plasma_complex32_t *C, int ldc);

void core_cher2k(plasma_enum_t uplo, plasma_enum_t trans,
                 int n, int k,
                 plasma_complex32_t alpha, const plasma_complex32_t *A, int lda,
                                           const plasma_complex32_t *B, int ldb,
                 float beta,                    plasma_complex32_t *C, int ldc);

void core_cherk(plasma_enum_t uplo, plasma_enum_t trans,
                int n, int k,
                float alpha, const plasma_complex32_t *A, int lda,
                float beta,        plasma_complex32_t *C, int ldc);

void core_chessq(plasma_enum_t uplo,
                 int n,
                 const plasma_complex32_t *A, int lda,
                 float *scale, float *sumsq);

void core_csyssq(plasma_enum_t uplo,
                 int n,
                 const plasma_complex32_t *A, int lda,
                 float *scale, float *sumsq);

void core_clacpy(plasma_enum_t uplo,
                 int m, int n,
                 const plasma_complex32_t *A, int lda,
                       plasma_complex32_t *B, int ldb);

void core_clacpy_lapack2tile_band(plasma_enum_t uplo,
                                  int it, int jt,
                                  int m, int n, int nb, int kl, int ku,
                                  const plasma_complex32_t *A, int lda,
                                        plasma_complex32_t *B, int ldb);

void core_clacpy_tile2lapack_band(plasma_enum_t uplo,
                                  int it, int jt,
                                  int m, int n, int nb, int kl, int ku,
                                  const plasma_complex32_t *B, int ldb,
                                        plasma_complex32_t *A, int lda);

void core_clange(plasma_enum_t norm,
                 int m, int n,
                 const plasma_complex32_t *A, int lda,
                 float *work, float *result);

void core_clanhe(plasma_enum_t norm, plasma_enum_t uplo,
                 int n,
                 const plasma_complex32_t *A, int lda,
                 float *work, float *value);

void core_clansy(plasma_enum_t norm, plasma_enum_t uplo,
                 int n,
                 const plasma_complex32_t *A, int lda,
                 float *work, float *value);

void core_clantr(plasma_enum_t norm, plasma_enum_t uplo, plasma_enum_t diag,
                 int m, int n,
                 const plasma_complex32_t *A, int lda,
                 float *work, float *value);

void core_clascl(plasma_enum_t uplo,
                 float cfrom, float cto,
                 int m, int n,
                 plasma_complex32_t *A, int lda);

void core_claset(plasma_enum_t uplo,
                 int m, int n,
                 plasma_complex32_t alpha, plasma_complex32_t beta,
                 plasma_complex32_t *A, int lda);

void core_claswp(plasma_enum_t colrow,
                 plasma_desc_t A, int k1, int k2, const int *ipiv, int incx);

int core_clauum(plasma_enum_t uplo,
                int n,
                plasma_complex32_t *A, int lda);

int core_cpamm(int op, plasma_enum_t side, plasma_enum_t storev,
               int m, int n, int k, int l,
               const plasma_complex32_t *A1, int lda1,
                     plasma_complex32_t *A2, int lda2,
               const plasma_complex32_t *V,  int ldv,
                     plasma_complex32_t *W,  int ldw);

int core_cparfb(plasma_enum_t side, plasma_enum_t trans, plasma_enum_t direct,
                plasma_enum_t storev,
                int m1, int n1, int m2, int n2, int k, int l,
                      plasma_complex32_t *A1,   int lda1,
                      plasma_complex32_t *A2,   int lda2,
                const plasma_complex32_t *V,    int ldv,
                const plasma_complex32_t *T,    int ldt,
                      plasma_complex32_t *work, int ldwork);

int core_cpemv(plasma_enum_t trans, int storev,
               int m, int n, int l,
               plasma_complex32_t alpha,
               const plasma_complex32_t *A, int lda,
               const plasma_complex32_t *X, int incx,
               plasma_complex32_t beta,
               plasma_complex32_t *Y, int incy,
               plasma_complex32_t *work);

int core_cpotrf(plasma_enum_t uplo,
                int n,
                plasma_complex32_t *A, int lda);

void core_csymm(plasma_enum_t side, plasma_enum_t uplo,
                int m, int n,
                plasma_complex32_t alpha, const plasma_complex32_t *A, int lda,
                                          const plasma_complex32_t *B, int ldb,
                plasma_complex32_t beta,        plasma_complex32_t *C, int ldc);

void core_csyr2k(
    plasma_enum_t uplo, plasma_enum_t trans,
    int n, int k,
    plasma_complex32_t alpha, const plasma_complex32_t *A, int lda,
                              const plasma_complex32_t *B, int ldb,
    plasma_complex32_t beta,        plasma_complex32_t *C, int ldc);

void core_csyrk(plasma_enum_t uplo, plasma_enum_t trans,
                int n, int k,
                plasma_complex32_t alpha, const plasma_complex32_t *A, int lda,
                plasma_complex32_t beta,        plasma_complex32_t *C, int ldc);

int core_ctradd(plasma_enum_t uplo, plasma_enum_t transa,
                int m, int n,
                plasma_complex32_t alpha, const plasma_complex32_t *A, int lda,
                plasma_complex32_t beta,        plasma_complex32_t *B, int ldb);

void core_ctrmm(plasma_enum_t side, plasma_enum_t uplo,
                plasma_enum_t transa, plasma_enum_t diag,
                int m, int n,
                plasma_complex32_t alpha, const plasma_complex32_t *A, int lda,
                                                plasma_complex32_t *B, int ldb);

void core_ctrsm(plasma_enum_t side, plasma_enum_t uplo,
                plasma_enum_t transa, plasma_enum_t diag,
                int m, int n,
                plasma_complex32_t alpha, const plasma_complex32_t *A, int lda,
                                                plasma_complex32_t *B, int ldb);

void core_ctrssq(plasma_enum_t uplo, plasma_enum_t diag,
                 int m, int n,
                 const plasma_complex32_t *A, int lda,
                 float *scale, float *sumsq);

int core_ctrtri(plasma_enum_t uplo, plasma_enum_t diag,
                int n,
                plasma_complex32_t *A, int lda);

int core_ctslqt(int m, int n, int ib,
                plasma_complex32_t *A1, int lda1,
                plasma_complex32_t *A2, int lda2,
                plasma_complex32_t *T,  int ldt,
                plasma_complex32_t *tau,
                plasma_complex32_t *work);

int core_ctsmlq(plasma_enum_t side, plasma_enum_t trans,
                int m1, int n1, int m2, int n2, int k, int ib,
                      plasma_complex32_t *A1,   int lda1,
                      plasma_complex32_t *A2,   int lda2,
                const plasma_complex32_t *V,    int ldv,
                const plasma_complex32_t *T,    int ldt,
                      plasma_complex32_t *work, int ldwork);

int core_ctsmqr(plasma_enum_t side, plasma_enum_t trans,
                int m1, int n1, int m2, int n2, int k, int ib,
                      plasma_complex32_t *A1,   int lda1,
                      plasma_complex32_t *A2,   int lda2,
                const plasma_complex32_t *V,    int ldv,
                const plasma_complex32_t *T,    int ldt,
                      plasma_complex32_t *work, int ldwork);

int core_ctsqrt(int m, int n, int ib,
                plasma_complex32_t *A1, int lda1,
                plasma_complex32_t *A2, int lda2,
                plasma_complex32_t *T,  int ldt,
                plasma_complex32_t *tau,
                plasma_complex32_t *work);

int core_cttlqt(int m, int n, int ib,
                plasma_complex32_t *A1, int lda1,
                plasma_complex32_t *A2, int lda2,
                plasma_complex32_t *T,  int ldt,
                plasma_complex32_t *tau,
                plasma_complex32_t *work);

int core_cttmlq(plasma_enum_t side, plasma_enum_t trans,
                int m1, int n1, int m2, int n2, int k, int ib,
                      plasma_complex32_t *A1,   int lda1,
                      plasma_complex32_t *A2,   int lda2,
                const plasma_complex32_t *V,    int ldv,
                const plasma_complex32_t *T,    int ldt,
                      plasma_complex32_t *work, int ldwork);

int core_cttmqr(plasma_enum_t side, plasma_enum_t trans,
                int m1, int n1, int m2, int n2, int k, int ib,
                      plasma_complex32_t *A1,   int lda1,
                      plasma_complex32_t *A2,   int lda2,
                const plasma_complex32_t *V,    int ldv,
                const plasma_complex32_t *T,    int ldt,
                      plasma_complex32_t *work, int ldwork);

int core_cttqrt(int m, int n, int ib,
                plasma_complex32_t *A1, int lda1,
                plasma_complex32_t *A2, int lda2,
                plasma_complex32_t *T,  int ldt,
                plasma_complex32_t *tau,
                plasma_complex32_t *work);

int core_cunmlq(plasma_enum_t side, plasma_enum_t trans,
                int m, int n, int k, int ib,
                const plasma_complex32_t *A,    int lda,
                const plasma_complex32_t *T,    int ldt,
                      plasma_complex32_t *C,    int ldc,
                      plasma_complex32_t *work, int ldwork);

int core_cunmqr(plasma_enum_t side, plasma_enum_t trans,
                int m, int n, int k, int ib,
                const plasma_complex32_t *A,    int lda,
                const plasma_complex32_t *T,    int ldt,
                      plasma_complex32_t *C,    int ldc,
                      plasma_complex32_t *work, int ldwork);

/******************************************************************************/
void core_omp_scamax(int colrow, int m, int n,
                     const plasma_complex32_t *A, int lda,
                     float *values,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_cgeadd(
    plasma_enum_t transa, int m, int n,
    plasma_complex32_t alpha, const plasma_complex32_t *A, int lda,
    plasma_complex32_t beta,        plasma_complex32_t *B, int ldb,
    plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_cgelqt(int m, int n, int ib,
                     plasma_complex32_t *A, int lda,
                     plasma_complex32_t *T, int ldt,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_cgemm(
    plasma_enum_t transa, plasma_enum_t transb,
    int m, int n, int k,
    plasma_complex32_t alpha, const plasma_complex32_t *A, int lda,
                              const plasma_complex32_t *B, int ldb,
    plasma_complex32_t beta,        plasma_complex32_t *C, int ldc,
    plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_cgeqrt(int m, int n, int ib,
                     plasma_complex32_t *A, int lda,
                     plasma_complex32_t *T, int ldt,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_cgessq(int m, int n,
                     const plasma_complex32_t *A, int lda,
                     float *scale, float *sumsq,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_cgessq_aux(int n,
                         const float *scale, const float *sumsq,
                         float *value,
                         plasma_sequence_t *sequence,
                         plasma_request_t *request);

void core_omp_chemm(
    plasma_enum_t side, plasma_enum_t uplo,
    int m, int n,
    plasma_complex32_t alpha, const plasma_complex32_t *A, int lda,
                              const plasma_complex32_t *B, int ldb,
    plasma_complex32_t beta,        plasma_complex32_t *C, int ldc,
    plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_cher2k(
    plasma_enum_t uplo, plasma_enum_t trans,
    int n, int k,
    plasma_complex32_t alpha, const plasma_complex32_t *A, int lda,
                              const plasma_complex32_t *B, int ldb,
    float beta,                    plasma_complex32_t *C, int ldc,
    plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_cherk(plasma_enum_t uplo, plasma_enum_t trans,
                    int n, int k,
                    float alpha, const plasma_complex32_t *A, int lda,
                    float beta,        plasma_complex32_t *C, int ldc,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_chessq(plasma_enum_t uplo,
                     int n,
                     const plasma_complex32_t *A, int lda,
                     float *scale, float *sumsq,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_csyssq(plasma_enum_t uplo,
                     int n,
                     const plasma_complex32_t *A, int lda,
                     float *scale, float *sumsq,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_csyssq_aux(int m, int n,
                         const float *scale, const float *sumsq,
                         float *value,
                         plasma_sequence_t *sequence,
                         plasma_request_t *request);

void core_omp_clacpy(plasma_enum_t uplo,
                     int m, int n,
                     const plasma_complex32_t *A, int lda,
                           plasma_complex32_t *B, int ldb,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_clacpy_lapack2tile_band(plasma_enum_t uplo,
                                      int it, int jt,
                                      int m, int n, int nb, int kl, int ku,
                                      const plasma_complex32_t *A, int lda,
                                            plasma_complex32_t *B, int ldb);

void core_omp_clacpy_tile2lapack_band(plasma_enum_t uplo,
                                      int it, int jt,
                                      int m, int n, int nb, int kl, int ku,
                                      const plasma_complex32_t *B, int ldb,
                                            plasma_complex32_t *A, int lda);

void core_omp_clange(plasma_enum_t norm,
                     int m, int n,
                     const plasma_complex32_t *A, int lda,
                     float *work, float *result,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_clange_aux(plasma_enum_t norm,
                         int m, int n,
                         const plasma_complex32_t *A, int lda,
                         float *value,
                         plasma_sequence_t *sequence,
                         plasma_request_t *request);

void core_omp_clanhe(plasma_enum_t norm, plasma_enum_t uplo,
                     int n,
                     const plasma_complex32_t *A, int lda,
                     float *work, float *value,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_clanhe_aux(plasma_enum_t norm, plasma_enum_t uplo,
                         int n,
                         const plasma_complex32_t *A, int lda,
                         float *value,
                         plasma_sequence_t *sequence,
                         plasma_request_t *request);

void core_omp_clansy(plasma_enum_t norm, plasma_enum_t uplo,
                     int n,
                     const plasma_complex32_t *A, int lda,
                     float *work, float *value,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_clansy_aux(plasma_enum_t norm, plasma_enum_t uplo,
                         int n,
                         const plasma_complex32_t *A, int lda,
                         float *value,
                         plasma_sequence_t *sequence,
                         plasma_request_t *request);

void core_omp_clantr(plasma_enum_t norm, plasma_enum_t uplo, plasma_enum_t diag,
                     int m, int n,
                     const plasma_complex32_t *A, int lda,
                     float *work, float *value,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_clantr_aux(plasma_enum_t norm, plasma_enum_t uplo,
                         plasma_enum_t diag,
                         int m, int n,
                         const plasma_complex32_t *A, int lda,
                         float *value,
                         plasma_sequence_t *sequence,
                         plasma_request_t *request);

void core_omp_clascl(plasma_enum_t uplo,
                     float cfrom, float cto,
                     int m, int n,
                     plasma_complex32_t *A, int lda,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_claset(plasma_enum_t uplo,
                     int mb, int nb,
                     int i, int j,
                     int m, int n,
                     plasma_complex32_t alpha, plasma_complex32_t beta,
                     plasma_complex32_t *A);

void core_omp_clauum(plasma_enum_t uplo,
                     int n,
                     plasma_complex32_t *A, int lda,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_cpotrf(plasma_enum_t uplo,
                     int n,
                     plasma_complex32_t *A, int lda,
                     int iinfo,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_csymm(
    plasma_enum_t side, plasma_enum_t uplo,
    int m, int n,
    plasma_complex32_t alpha, const plasma_complex32_t *A, int lda,
                              const plasma_complex32_t *B, int ldb,
    plasma_complex32_t beta,        plasma_complex32_t *C, int ldc,
    plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_csyr2k(
    plasma_enum_t uplo, plasma_enum_t trans,
    int n, int k,
    plasma_complex32_t alpha, const plasma_complex32_t *A, int lda,
                              const plasma_complex32_t *B, int ldb,
    plasma_complex32_t beta,        plasma_complex32_t *C, int ldc,
    plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_csyrk(
    plasma_enum_t uplo, plasma_enum_t trans,
    int n, int k,
    plasma_complex32_t alpha, const plasma_complex32_t *A, int lda,
    plasma_complex32_t beta,        plasma_complex32_t *C, int ldc,
    plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_ctradd(
    plasma_enum_t uplo, plasma_enum_t transa,
    int m, int n,
    plasma_complex32_t alpha, const plasma_complex32_t *A, int lda,
    plasma_complex32_t beta,        plasma_complex32_t *B, int ldb,
    plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_ctrmm(
    plasma_enum_t side, plasma_enum_t uplo,
    plasma_enum_t transa, plasma_enum_t diag,
    int m, int n,
    plasma_complex32_t alpha, const plasma_complex32_t *A, int lda,
                                    plasma_complex32_t *B, int ldb,
    plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_ctrsm(
    plasma_enum_t side, plasma_enum_t uplo,
    plasma_enum_t transa, plasma_enum_t diag,
    int m, int n,
    plasma_complex32_t alpha, const plasma_complex32_t *A, int lda,
                                    plasma_complex32_t *B, int ldb,
    plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_ctrssq(plasma_enum_t uplo, plasma_enum_t diag,
                     int m, int n,
                     const plasma_complex32_t *A, int lda,
                     float *scale, float *sumsq,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_ctrtri(plasma_enum_t uplo, plasma_enum_t diag,
                     int n,
                     plasma_complex32_t *A, int lda,
                     int iinfo,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_ctslqt(int m, int n, int ib,
                     plasma_complex32_t *A1, int lda1,
                     plasma_complex32_t *A2, int lda2,
                     plasma_complex32_t *T,  int ldt,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_ctsmlq(plasma_enum_t side, plasma_enum_t trans,
                     int m1, int n1, int m2, int n2, int k, int ib,
                           plasma_complex32_t *A1, int lda1,
                           plasma_complex32_t *A2, int lda2,
                     const plasma_complex32_t *V,  int ldv,
                     const plasma_complex32_t *T,  int ldt,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_ctsmqr(plasma_enum_t side, plasma_enum_t trans,
                     int m1, int n1, int m2, int n2, int k, int ib,
                           plasma_complex32_t *A1, int lda1,
                           plasma_complex32_t *A2, int lda2,
                     const plasma_complex32_t *V, int ldv,
                     const plasma_complex32_t *T, int ldt,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_ctsqrt(int m, int n, int ib,
                     plasma_complex32_t *A1, int lda1,
                     plasma_complex32_t *A2, int lda2,
                     plasma_complex32_t *T,  int ldt,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_cttlqt(int m, int n, int ib,
                     plasma_complex32_t *A1, int lda1,
                     plasma_complex32_t *A2, int lda2,
                     plasma_complex32_t *T,  int ldt,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_cttmlq(plasma_enum_t side, plasma_enum_t trans,
                     int m1, int n1, int m2, int n2, int k, int ib,
                           plasma_complex32_t *A1, int lda1,
                           plasma_complex32_t *A2, int lda2,
                     const plasma_complex32_t *V,  int ldv,
                     const plasma_complex32_t *T,  int ldt,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_cttmqr(plasma_enum_t side, plasma_enum_t trans,
                     int m1, int n1, int m2, int n2, int k, int ib,
                           plasma_complex32_t *A1, int lda1,
                           plasma_complex32_t *A2, int lda2,
                     const plasma_complex32_t *V, int ldv,
                     const plasma_complex32_t *T, int ldt,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_cttqrt(int m, int n, int ib,
                     plasma_complex32_t *A1, int lda1,
                     plasma_complex32_t *A2, int lda2,
                     plasma_complex32_t *T,  int ldt,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_cunmlq(plasma_enum_t side, plasma_enum_t trans,
                     int m, int n, int k, int ib,
                     const plasma_complex32_t *A, int lda,
                     const plasma_complex32_t *T, int ldt,
                           plasma_complex32_t *C, int ldc,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_cunmqr(plasma_enum_t side, plasma_enum_t trans,
                     int m, int n, int k, int ib,
                     const plasma_complex32_t *A, int lda,
                     const plasma_complex32_t *T, int ldt,
                           plasma_complex32_t *C, int ldc,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

#undef COMPLEX

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // ICL_CORE_BLAS_C_H
