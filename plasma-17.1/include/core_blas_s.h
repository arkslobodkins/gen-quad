/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from include/core_blas_z.h, normal z -> s, Mon Nov 22 19:09:33 2021
 *
 **/
#ifndef ICL_CORE_BLAS_S_H
#define ICL_CORE_BLAS_S_H

#include "plasma_async.h"
#include "plasma_barrier.h"
#include "plasma_descriptor.h"
#include "plasma_types.h"
#include "plasma_workspace.h"
#include "plasma_descriptor.h"

#ifdef __cplusplus
extern "C" {
#endif

#define REAL

/******************************************************************************/
#ifdef COMPLEX
float fabsf(float alpha);
#endif

int core_sgeadd(plasma_enum_t transa,
                int m, int n,
                float alpha, const float *A, int lda,
                float beta,        float *B, int ldb);

int core_sgelqt(int m, int n, int ib,
                float *A, int lda,
                float *T, int ldt,
                float *tau,
                float *work);

void core_sgemm(plasma_enum_t transa, plasma_enum_t transb,
                int m, int n, int k,
                float alpha, const float *A, int lda,
                                          const float *B, int ldb,
                float beta,        float *C, int ldc);

int core_sgeqrt(int m, int n, int ib,
                float *A, int lda,
                float *T, int ldt,
                float *tau,
                float *work);

void core_sgessq(int m, int n,
                 const float *A, int lda,
                 float *scale, float *sumsq);

int core_sgetrf(plasma_desc_t A, int *ipiv, int ib, int rank, int size,
                plasma_barrier_t *barrier);

void core_ssymm(plasma_enum_t side, plasma_enum_t uplo,
                int m, int n,
                float alpha, const float *A, int lda,
                                          const float *B, int ldb,
                float beta,        float *C, int ldc);

void core_ssyr2k(plasma_enum_t uplo, plasma_enum_t trans,
                 int n, int k,
                 float alpha, const float *A, int lda,
                                           const float *B, int ldb,
                 float beta,                    float *C, int ldc);

void core_ssyrk(plasma_enum_t uplo, plasma_enum_t trans,
                int n, int k,
                float alpha, const float *A, int lda,
                float beta,        float *C, int ldc);

void core_ssyssq(plasma_enum_t uplo,
                 int n,
                 const float *A, int lda,
                 float *scale, float *sumsq);

void core_ssyssq(plasma_enum_t uplo,
                 int n,
                 const float *A, int lda,
                 float *scale, float *sumsq);

void core_slacpy(plasma_enum_t uplo,
                 int m, int n,
                 const float *A, int lda,
                       float *B, int ldb);

void core_slacpy_lapack2tile_band(plasma_enum_t uplo,
                                  int it, int jt,
                                  int m, int n, int nb, int kl, int ku,
                                  const float *A, int lda,
                                        float *B, int ldb);

void core_slacpy_tile2lapack_band(plasma_enum_t uplo,
                                  int it, int jt,
                                  int m, int n, int nb, int kl, int ku,
                                  const float *B, int ldb,
                                        float *A, int lda);

void core_slange(plasma_enum_t norm,
                 int m, int n,
                 const float *A, int lda,
                 float *work, float *result);

void core_slansy(plasma_enum_t norm, plasma_enum_t uplo,
                 int n,
                 const float *A, int lda,
                 float *work, float *value);

void core_slansy(plasma_enum_t norm, plasma_enum_t uplo,
                 int n,
                 const float *A, int lda,
                 float *work, float *value);

void core_slantr(plasma_enum_t norm, plasma_enum_t uplo, plasma_enum_t diag,
                 int m, int n,
                 const float *A, int lda,
                 float *work, float *value);

void core_slascl(plasma_enum_t uplo,
                 float cfrom, float cto,
                 int m, int n,
                 float *A, int lda);

void core_slaset(plasma_enum_t uplo,
                 int m, int n,
                 float alpha, float beta,
                 float *A, int lda);

void core_slaswp(plasma_enum_t colrow,
                 plasma_desc_t A, int k1, int k2, const int *ipiv, int incx);

int core_slauum(plasma_enum_t uplo,
                int n,
                float *A, int lda);

int core_spamm(int op, plasma_enum_t side, plasma_enum_t storev,
               int m, int n, int k, int l,
               const float *A1, int lda1,
                     float *A2, int lda2,
               const float *V,  int ldv,
                     float *W,  int ldw);

int core_sparfb(plasma_enum_t side, plasma_enum_t trans, plasma_enum_t direct,
                plasma_enum_t storev,
                int m1, int n1, int m2, int n2, int k, int l,
                      float *A1,   int lda1,
                      float *A2,   int lda2,
                const float *V,    int ldv,
                const float *T,    int ldt,
                      float *work, int ldwork);

int core_spemv(plasma_enum_t trans, int storev,
               int m, int n, int l,
               float alpha,
               const float *A, int lda,
               const float *X, int incx,
               float beta,
               float *Y, int incy,
               float *work);

int core_spotrf(plasma_enum_t uplo,
                int n,
                float *A, int lda);

void core_ssymm(plasma_enum_t side, plasma_enum_t uplo,
                int m, int n,
                float alpha, const float *A, int lda,
                                          const float *B, int ldb,
                float beta,        float *C, int ldc);

void core_ssyr2k(
    plasma_enum_t uplo, plasma_enum_t trans,
    int n, int k,
    float alpha, const float *A, int lda,
                              const float *B, int ldb,
    float beta,        float *C, int ldc);

void core_ssyrk(plasma_enum_t uplo, plasma_enum_t trans,
                int n, int k,
                float alpha, const float *A, int lda,
                float beta,        float *C, int ldc);

int core_stradd(plasma_enum_t uplo, plasma_enum_t transa,
                int m, int n,
                float alpha, const float *A, int lda,
                float beta,        float *B, int ldb);

void core_strmm(plasma_enum_t side, plasma_enum_t uplo,
                plasma_enum_t transa, plasma_enum_t diag,
                int m, int n,
                float alpha, const float *A, int lda,
                                                float *B, int ldb);

void core_strsm(plasma_enum_t side, plasma_enum_t uplo,
                plasma_enum_t transa, plasma_enum_t diag,
                int m, int n,
                float alpha, const float *A, int lda,
                                                float *B, int ldb);

void core_strssq(plasma_enum_t uplo, plasma_enum_t diag,
                 int m, int n,
                 const float *A, int lda,
                 float *scale, float *sumsq);

int core_strtri(plasma_enum_t uplo, plasma_enum_t diag,
                int n,
                float *A, int lda);

int core_stslqt(int m, int n, int ib,
                float *A1, int lda1,
                float *A2, int lda2,
                float *T,  int ldt,
                float *tau,
                float *work);

int core_stsmlq(plasma_enum_t side, plasma_enum_t trans,
                int m1, int n1, int m2, int n2, int k, int ib,
                      float *A1,   int lda1,
                      float *A2,   int lda2,
                const float *V,    int ldv,
                const float *T,    int ldt,
                      float *work, int ldwork);

int core_stsmqr(plasma_enum_t side, plasma_enum_t trans,
                int m1, int n1, int m2, int n2, int k, int ib,
                      float *A1,   int lda1,
                      float *A2,   int lda2,
                const float *V,    int ldv,
                const float *T,    int ldt,
                      float *work, int ldwork);

int core_stsqrt(int m, int n, int ib,
                float *A1, int lda1,
                float *A2, int lda2,
                float *T,  int ldt,
                float *tau,
                float *work);

int core_sttlqt(int m, int n, int ib,
                float *A1, int lda1,
                float *A2, int lda2,
                float *T,  int ldt,
                float *tau,
                float *work);

int core_sttmlq(plasma_enum_t side, plasma_enum_t trans,
                int m1, int n1, int m2, int n2, int k, int ib,
                      float *A1,   int lda1,
                      float *A2,   int lda2,
                const float *V,    int ldv,
                const float *T,    int ldt,
                      float *work, int ldwork);

int core_sttmqr(plasma_enum_t side, plasma_enum_t trans,
                int m1, int n1, int m2, int n2, int k, int ib,
                      float *A1,   int lda1,
                      float *A2,   int lda2,
                const float *V,    int ldv,
                const float *T,    int ldt,
                      float *work, int ldwork);

int core_sttqrt(int m, int n, int ib,
                float *A1, int lda1,
                float *A2, int lda2,
                float *T,  int ldt,
                float *tau,
                float *work);

int core_sormlq(plasma_enum_t side, plasma_enum_t trans,
                int m, int n, int k, int ib,
                const float *A,    int lda,
                const float *T,    int ldt,
                      float *C,    int ldc,
                      float *work, int ldwork);

int core_sormqr(plasma_enum_t side, plasma_enum_t trans,
                int m, int n, int k, int ib,
                const float *A,    int lda,
                const float *T,    int ldt,
                      float *C,    int ldc,
                      float *work, int ldwork);

/******************************************************************************/
void core_omp_samax(int colrow, int m, int n,
                     const float *A, int lda,
                     float *values,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_sgeadd(
    plasma_enum_t transa, int m, int n,
    float alpha, const float *A, int lda,
    float beta,        float *B, int ldb,
    plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_sgelqt(int m, int n, int ib,
                     float *A, int lda,
                     float *T, int ldt,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_sgemm(
    plasma_enum_t transa, plasma_enum_t transb,
    int m, int n, int k,
    float alpha, const float *A, int lda,
                              const float *B, int ldb,
    float beta,        float *C, int ldc,
    plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_sgeqrt(int m, int n, int ib,
                     float *A, int lda,
                     float *T, int ldt,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_sgessq(int m, int n,
                     const float *A, int lda,
                     float *scale, float *sumsq,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_sgessq_aux(int n,
                         const float *scale, const float *sumsq,
                         float *value,
                         plasma_sequence_t *sequence,
                         plasma_request_t *request);

void core_omp_ssymm(
    plasma_enum_t side, plasma_enum_t uplo,
    int m, int n,
    float alpha, const float *A, int lda,
                              const float *B, int ldb,
    float beta,        float *C, int ldc,
    plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_ssyr2k(
    plasma_enum_t uplo, plasma_enum_t trans,
    int n, int k,
    float alpha, const float *A, int lda,
                              const float *B, int ldb,
    float beta,                    float *C, int ldc,
    plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_ssyrk(plasma_enum_t uplo, plasma_enum_t trans,
                    int n, int k,
                    float alpha, const float *A, int lda,
                    float beta,        float *C, int ldc,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_ssyssq(plasma_enum_t uplo,
                     int n,
                     const float *A, int lda,
                     float *scale, float *sumsq,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_ssyssq(plasma_enum_t uplo,
                     int n,
                     const float *A, int lda,
                     float *scale, float *sumsq,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_ssyssq_aux(int m, int n,
                         const float *scale, const float *sumsq,
                         float *value,
                         plasma_sequence_t *sequence,
                         plasma_request_t *request);

void core_omp_slacpy(plasma_enum_t uplo,
                     int m, int n,
                     const float *A, int lda,
                           float *B, int ldb,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_slacpy_lapack2tile_band(plasma_enum_t uplo,
                                      int it, int jt,
                                      int m, int n, int nb, int kl, int ku,
                                      const float *A, int lda,
                                            float *B, int ldb);

void core_omp_slacpy_tile2lapack_band(plasma_enum_t uplo,
                                      int it, int jt,
                                      int m, int n, int nb, int kl, int ku,
                                      const float *B, int ldb,
                                            float *A, int lda);

void core_omp_slange(plasma_enum_t norm,
                     int m, int n,
                     const float *A, int lda,
                     float *work, float *result,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_slange_aux(plasma_enum_t norm,
                         int m, int n,
                         const float *A, int lda,
                         float *value,
                         plasma_sequence_t *sequence,
                         plasma_request_t *request);

void core_omp_slansy(plasma_enum_t norm, plasma_enum_t uplo,
                     int n,
                     const float *A, int lda,
                     float *work, float *value,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_slansy_aux(plasma_enum_t norm, plasma_enum_t uplo,
                         int n,
                         const float *A, int lda,
                         float *value,
                         plasma_sequence_t *sequence,
                         plasma_request_t *request);

void core_omp_slansy(plasma_enum_t norm, plasma_enum_t uplo,
                     int n,
                     const float *A, int lda,
                     float *work, float *value,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_slansy_aux(plasma_enum_t norm, plasma_enum_t uplo,
                         int n,
                         const float *A, int lda,
                         float *value,
                         plasma_sequence_t *sequence,
                         plasma_request_t *request);

void core_omp_slantr(plasma_enum_t norm, plasma_enum_t uplo, plasma_enum_t diag,
                     int m, int n,
                     const float *A, int lda,
                     float *work, float *value,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_slantr_aux(plasma_enum_t norm, plasma_enum_t uplo,
                         plasma_enum_t diag,
                         int m, int n,
                         const float *A, int lda,
                         float *value,
                         plasma_sequence_t *sequence,
                         plasma_request_t *request);

void core_omp_slascl(plasma_enum_t uplo,
                     float cfrom, float cto,
                     int m, int n,
                     float *A, int lda,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_slaset(plasma_enum_t uplo,
                     int mb, int nb,
                     int i, int j,
                     int m, int n,
                     float alpha, float beta,
                     float *A);

void core_omp_slauum(plasma_enum_t uplo,
                     int n,
                     float *A, int lda,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_spotrf(plasma_enum_t uplo,
                     int n,
                     float *A, int lda,
                     int iinfo,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_ssymm(
    plasma_enum_t side, plasma_enum_t uplo,
    int m, int n,
    float alpha, const float *A, int lda,
                              const float *B, int ldb,
    float beta,        float *C, int ldc,
    plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_ssyr2k(
    plasma_enum_t uplo, plasma_enum_t trans,
    int n, int k,
    float alpha, const float *A, int lda,
                              const float *B, int ldb,
    float beta,        float *C, int ldc,
    plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_ssyrk(
    plasma_enum_t uplo, plasma_enum_t trans,
    int n, int k,
    float alpha, const float *A, int lda,
    float beta,        float *C, int ldc,
    plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_stradd(
    plasma_enum_t uplo, plasma_enum_t transa,
    int m, int n,
    float alpha, const float *A, int lda,
    float beta,        float *B, int ldb,
    plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_strmm(
    plasma_enum_t side, plasma_enum_t uplo,
    plasma_enum_t transa, plasma_enum_t diag,
    int m, int n,
    float alpha, const float *A, int lda,
                                    float *B, int ldb,
    plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_strsm(
    plasma_enum_t side, plasma_enum_t uplo,
    plasma_enum_t transa, plasma_enum_t diag,
    int m, int n,
    float alpha, const float *A, int lda,
                                    float *B, int ldb,
    plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_strssq(plasma_enum_t uplo, plasma_enum_t diag,
                     int m, int n,
                     const float *A, int lda,
                     float *scale, float *sumsq,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_strtri(plasma_enum_t uplo, plasma_enum_t diag,
                     int n,
                     float *A, int lda,
                     int iinfo,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_stslqt(int m, int n, int ib,
                     float *A1, int lda1,
                     float *A2, int lda2,
                     float *T,  int ldt,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_stsmlq(plasma_enum_t side, plasma_enum_t trans,
                     int m1, int n1, int m2, int n2, int k, int ib,
                           float *A1, int lda1,
                           float *A2, int lda2,
                     const float *V,  int ldv,
                     const float *T,  int ldt,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_stsmqr(plasma_enum_t side, plasma_enum_t trans,
                     int m1, int n1, int m2, int n2, int k, int ib,
                           float *A1, int lda1,
                           float *A2, int lda2,
                     const float *V, int ldv,
                     const float *T, int ldt,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_stsqrt(int m, int n, int ib,
                     float *A1, int lda1,
                     float *A2, int lda2,
                     float *T,  int ldt,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_sttlqt(int m, int n, int ib,
                     float *A1, int lda1,
                     float *A2, int lda2,
                     float *T,  int ldt,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_sttmlq(plasma_enum_t side, plasma_enum_t trans,
                     int m1, int n1, int m2, int n2, int k, int ib,
                           float *A1, int lda1,
                           float *A2, int lda2,
                     const float *V,  int ldv,
                     const float *T,  int ldt,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_sttmqr(plasma_enum_t side, plasma_enum_t trans,
                     int m1, int n1, int m2, int n2, int k, int ib,
                           float *A1, int lda1,
                           float *A2, int lda2,
                     const float *V, int ldv,
                     const float *T, int ldt,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_sttqrt(int m, int n, int ib,
                     float *A1, int lda1,
                     float *A2, int lda2,
                     float *T,  int ldt,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_sormlq(plasma_enum_t side, plasma_enum_t trans,
                     int m, int n, int k, int ib,
                     const float *A, int lda,
                     const float *T, int ldt,
                           float *C, int ldc,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

void core_omp_sormqr(plasma_enum_t side, plasma_enum_t trans,
                     int m, int n, int k, int ib,
                     const float *A, int lda,
                     const float *T, int ldt,
                           float *C, int ldc,
                     plasma_workspace_t work,
                     plasma_sequence_t *sequence, plasma_request_t *request);

#undef REAL

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // ICL_CORE_BLAS_S_H
