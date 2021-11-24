/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from include/plasma_internal_z.h, normal z -> c, Mon Nov 22 19:09:33 2021
 *
 **/
#ifndef ICL_PLASMA_INTERNAL_C_H
#define ICL_PLASMA_INTERNAL_C_H

#include "plasma_async.h"
#include "plasma_descriptor.h"
#include "plasma_types.h"
#include "plasma_workspace.h"

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************/
void plasma_pscamax(plasma_enum_t colrow,
                    plasma_desc_t A, float *work, float *values,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pcgbtrf(plasma_desc_t A, int *ipiv,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pcdesc2ge(plasma_desc_t A,
                      plasma_complex32_t *pA, int lda,
                      plasma_sequence_t *sequence,
                      plasma_request_t *request);

void plasma_pcdesc2pb(plasma_desc_t A,
                      plasma_complex32_t *pA, int lda,
                      plasma_sequence_t *sequence,
                      plasma_request_t *request);

void plasma_pcge2desc(plasma_complex32_t *pA, int lda,
                      plasma_desc_t A,
                      plasma_sequence_t *sequence,
                      plasma_request_t *request);

void plasma_pcgeadd(plasma_enum_t transa,
                    plasma_complex32_t alpha,  plasma_desc_t A,
                    plasma_complex32_t beta,   plasma_desc_t B,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pcgelqf(plasma_desc_t A, plasma_desc_t T,
                    plasma_workspace_t work,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pcgelqfrh(plasma_desc_t A, plasma_desc_t T,
                      plasma_workspace_t work,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pcgemm(plasma_enum_t transa, plasma_enum_t transb,
                   plasma_complex32_t alpha, plasma_desc_t A,
                                             plasma_desc_t B,
                   plasma_complex32_t beta,  plasma_desc_t C,
                   plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pcgeqrf(plasma_desc_t A, plasma_desc_t T,
                    plasma_workspace_t work,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pcgeqrfrh(plasma_desc_t A, plasma_desc_t T,
                      plasma_workspace_t work,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pcgetri_aux(plasma_desc_t A, plasma_desc_t W,
                        plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pcgetrf(plasma_desc_t A, int *ipiv,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pchemm(plasma_enum_t side, plasma_enum_t uplo,
                   plasma_complex32_t alpha, plasma_desc_t A,
                                             plasma_desc_t B,
                   plasma_complex32_t beta,  plasma_desc_t C,
                   plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pcher2k(plasma_enum_t uplo, plasma_enum_t trans,
                    plasma_complex32_t alpha, plasma_desc_t A,
                                              plasma_desc_t B,
                    float beta,              plasma_desc_t C,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pcherk(plasma_enum_t uplo, plasma_enum_t trans,
                   float alpha, plasma_desc_t A,
                   float beta,  plasma_desc_t C,
                   plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pclacpy(plasma_enum_t uplo, plasma_desc_t A, plasma_desc_t B,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pclange(plasma_enum_t norm,
                    plasma_desc_t A, float *work, float *value,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pclanhe(plasma_enum_t norm, plasma_enum_t uplo,
                    plasma_desc_t A, float *work, float *value,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pclansy(plasma_enum_t norm, plasma_enum_t uplo,
                    plasma_desc_t A, float *work, float *value,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pclantr(plasma_enum_t norm, plasma_enum_t uplo, plasma_enum_t diag,
                    plasma_desc_t A, float *work, float *value,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pclascl(plasma_enum_t uplo,
                    float cfrom, float cto,
                    plasma_desc_t A,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pclaset(plasma_enum_t uplo,
                    plasma_complex32_t alpha, plasma_complex32_t beta,
                    plasma_desc_t A,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pclaswp(plasma_enum_t colrow,
                    plasma_desc_t A, int *ipiv, int incx,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pclauum(plasma_enum_t uplo, plasma_desc_t A,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pcpb2desc(plasma_complex32_t *pA, int lda,
                      plasma_desc_t A,
                      plasma_sequence_t *sequence,
                      plasma_request_t *request);

void plasma_pcpbtrf(plasma_enum_t uplo, plasma_desc_t A,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pcpotrf(plasma_enum_t uplo, plasma_desc_t A,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pcsymm(plasma_enum_t side, plasma_enum_t uplo,
                   plasma_complex32_t alpha, plasma_desc_t A,
                                             plasma_desc_t B,
                   plasma_complex32_t beta,  plasma_desc_t C,
                   plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pcsyr2k(plasma_enum_t uplo, plasma_enum_t trans,
                    plasma_complex32_t alpha, plasma_desc_t A,
                                              plasma_desc_t B,
                    plasma_complex32_t beta,  plasma_desc_t C,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pcsyrk(plasma_enum_t uplo, plasma_enum_t trans,
                   plasma_complex32_t alpha, plasma_desc_t A,
                   plasma_complex32_t beta,  plasma_desc_t C,
                   plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pctbsm(plasma_enum_t side, plasma_enum_t uplo,
                   plasma_enum_t trans, plasma_enum_t diag,
                   plasma_complex32_t alpha, plasma_desc_t A,
                                             plasma_desc_t B,
                   const int *ipiv,
                   plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pctradd(plasma_enum_t uplo, plasma_enum_t transa,
                    plasma_complex32_t alpha,  plasma_desc_t A,
                    plasma_complex32_t beta,   plasma_desc_t B,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pctrmm(plasma_enum_t side, plasma_enum_t uplo,
                   plasma_enum_t trans, plasma_enum_t diag,
                   plasma_complex32_t alpha, plasma_desc_t A,
                                             plasma_desc_t B,
                   plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pctrsm(plasma_enum_t side, plasma_enum_t uplo,
                   plasma_enum_t trans, plasma_enum_t diag,
                   plasma_complex32_t alpha, plasma_desc_t A,
                                             plasma_desc_t B,
                   plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pctrtri(plasma_enum_t uplo, plasma_enum_t diag,
                    plasma_desc_t A,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pcunglq(plasma_desc_t A, plasma_desc_t T, plasma_desc_t Q,
                    plasma_workspace_t work,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pcunglqrh(plasma_desc_t A, plasma_desc_t T, plasma_desc_t Q,
                      plasma_workspace_t work,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pcungqr(plasma_desc_t A, plasma_desc_t T, plasma_desc_t Q,
                    plasma_workspace_t work,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pcungqrrh(plasma_desc_t A, plasma_desc_t T, plasma_desc_t Q,
                      plasma_workspace_t work,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pcunmlq(plasma_enum_t side, plasma_enum_t trans,
                    plasma_desc_t A, plasma_desc_t T, plasma_desc_t B,
                    plasma_workspace_t work,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pcunmlqrh(plasma_enum_t side, plasma_enum_t trans,
                      plasma_desc_t A, plasma_desc_t T, plasma_desc_t B,
                      plasma_workspace_t work,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pcunmqr(plasma_enum_t side, plasma_enum_t trans,
                    plasma_desc_t A, plasma_desc_t T, plasma_desc_t B,
                    plasma_workspace_t work,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pcunmqrrh(plasma_enum_t side, plasma_enum_t trans,
                      plasma_desc_t A, plasma_desc_t T, plasma_desc_t B,
                      plasma_workspace_t work,
                      plasma_sequence_t *sequence, plasma_request_t *request);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // ICL_PLASMA_INTERNAL_C_H
