/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from include/plasma_internal_z.h, normal z -> s, Mon Nov 22 19:09:33 2021
 *
 **/
#ifndef ICL_PLASMA_INTERNAL_S_H
#define ICL_PLASMA_INTERNAL_S_H

#include "plasma_async.h"
#include "plasma_descriptor.h"
#include "plasma_types.h"
#include "plasma_workspace.h"

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************/
void plasma_psamax(plasma_enum_t colrow,
                    plasma_desc_t A, float *work, float *values,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_psgbtrf(plasma_desc_t A, int *ipiv,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_psdesc2ge(plasma_desc_t A,
                      float *pA, int lda,
                      plasma_sequence_t *sequence,
                      plasma_request_t *request);

void plasma_psdesc2pb(plasma_desc_t A,
                      float *pA, int lda,
                      plasma_sequence_t *sequence,
                      plasma_request_t *request);

void plasma_psge2desc(float *pA, int lda,
                      plasma_desc_t A,
                      plasma_sequence_t *sequence,
                      plasma_request_t *request);

void plasma_psgeadd(plasma_enum_t transa,
                    float alpha,  plasma_desc_t A,
                    float beta,   plasma_desc_t B,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_psgelqf(plasma_desc_t A, plasma_desc_t T,
                    plasma_workspace_t work,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_psgelqfrh(plasma_desc_t A, plasma_desc_t T,
                      plasma_workspace_t work,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_psgemm(plasma_enum_t transa, plasma_enum_t transb,
                   float alpha, plasma_desc_t A,
                                             plasma_desc_t B,
                   float beta,  plasma_desc_t C,
                   plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_psgeqrf(plasma_desc_t A, plasma_desc_t T,
                    plasma_workspace_t work,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_psgeqrfrh(plasma_desc_t A, plasma_desc_t T,
                      plasma_workspace_t work,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_psgetri_aux(plasma_desc_t A, plasma_desc_t W,
                        plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_psgetrf(plasma_desc_t A, int *ipiv,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pssymm(plasma_enum_t side, plasma_enum_t uplo,
                   float alpha, plasma_desc_t A,
                                             plasma_desc_t B,
                   float beta,  plasma_desc_t C,
                   plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pssyr2k(plasma_enum_t uplo, plasma_enum_t trans,
                    float alpha, plasma_desc_t A,
                                              plasma_desc_t B,
                    float beta,              plasma_desc_t C,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pssyrk(plasma_enum_t uplo, plasma_enum_t trans,
                   float alpha, plasma_desc_t A,
                   float beta,  plasma_desc_t C,
                   plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pslacpy(plasma_enum_t uplo, plasma_desc_t A, plasma_desc_t B,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pslange(plasma_enum_t norm,
                    plasma_desc_t A, float *work, float *value,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pslansy(plasma_enum_t norm, plasma_enum_t uplo,
                    plasma_desc_t A, float *work, float *value,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pslansy(plasma_enum_t norm, plasma_enum_t uplo,
                    plasma_desc_t A, float *work, float *value,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pslantr(plasma_enum_t norm, plasma_enum_t uplo, plasma_enum_t diag,
                    plasma_desc_t A, float *work, float *value,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pslascl(plasma_enum_t uplo,
                    float cfrom, float cto,
                    plasma_desc_t A,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pslaset(plasma_enum_t uplo,
                    float alpha, float beta,
                    plasma_desc_t A,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pslaswp(plasma_enum_t colrow,
                    plasma_desc_t A, int *ipiv, int incx,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pslauum(plasma_enum_t uplo, plasma_desc_t A,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pspb2desc(float *pA, int lda,
                      plasma_desc_t A,
                      plasma_sequence_t *sequence,
                      plasma_request_t *request);

void plasma_pspbtrf(plasma_enum_t uplo, plasma_desc_t A,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pspotrf(plasma_enum_t uplo, plasma_desc_t A,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pssymm(plasma_enum_t side, plasma_enum_t uplo,
                   float alpha, plasma_desc_t A,
                                             plasma_desc_t B,
                   float beta,  plasma_desc_t C,
                   plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pssyr2k(plasma_enum_t uplo, plasma_enum_t trans,
                    float alpha, plasma_desc_t A,
                                              plasma_desc_t B,
                    float beta,  plasma_desc_t C,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pssyrk(plasma_enum_t uplo, plasma_enum_t trans,
                   float alpha, plasma_desc_t A,
                   float beta,  plasma_desc_t C,
                   plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pstbsm(plasma_enum_t side, plasma_enum_t uplo,
                   plasma_enum_t trans, plasma_enum_t diag,
                   float alpha, plasma_desc_t A,
                                             plasma_desc_t B,
                   const int *ipiv,
                   plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pstradd(plasma_enum_t uplo, plasma_enum_t transa,
                    float alpha,  plasma_desc_t A,
                    float beta,   plasma_desc_t B,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pstrmm(plasma_enum_t side, plasma_enum_t uplo,
                   plasma_enum_t trans, plasma_enum_t diag,
                   float alpha, plasma_desc_t A,
                                             plasma_desc_t B,
                   plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pstrsm(plasma_enum_t side, plasma_enum_t uplo,
                   plasma_enum_t trans, plasma_enum_t diag,
                   float alpha, plasma_desc_t A,
                                             plasma_desc_t B,
                   plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pstrtri(plasma_enum_t uplo, plasma_enum_t diag,
                    plasma_desc_t A,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_psorglq(plasma_desc_t A, plasma_desc_t T, plasma_desc_t Q,
                    plasma_workspace_t work,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_psorglqrh(plasma_desc_t A, plasma_desc_t T, plasma_desc_t Q,
                      plasma_workspace_t work,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_psorgqr(plasma_desc_t A, plasma_desc_t T, plasma_desc_t Q,
                    plasma_workspace_t work,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_psorgqrrh(plasma_desc_t A, plasma_desc_t T, plasma_desc_t Q,
                      plasma_workspace_t work,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_psormlq(plasma_enum_t side, plasma_enum_t trans,
                    plasma_desc_t A, plasma_desc_t T, plasma_desc_t B,
                    plasma_workspace_t work,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_psormlqrh(plasma_enum_t side, plasma_enum_t trans,
                      plasma_desc_t A, plasma_desc_t T, plasma_desc_t B,
                      plasma_workspace_t work,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_psormqr(plasma_enum_t side, plasma_enum_t trans,
                    plasma_desc_t A, plasma_desc_t T, plasma_desc_t B,
                    plasma_workspace_t work,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_psormqrrh(plasma_enum_t side, plasma_enum_t trans,
                      plasma_desc_t A, plasma_desc_t T, plasma_desc_t B,
                      plasma_workspace_t work,
                      plasma_sequence_t *sequence, plasma_request_t *request);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // ICL_PLASMA_INTERNAL_S_H
