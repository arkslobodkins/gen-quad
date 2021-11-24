/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from include/plasma_internal_z.h, normal z -> d, Mon Nov 22 19:09:33 2021
 *
 **/
#ifndef ICL_PLASMA_INTERNAL_D_H
#define ICL_PLASMA_INTERNAL_D_H

#include "plasma_async.h"
#include "plasma_descriptor.h"
#include "plasma_types.h"
#include "plasma_workspace.h"

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************/
void plasma_pdamax(plasma_enum_t colrow,
                    plasma_desc_t A, double *work, double *values,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pdgbtrf(plasma_desc_t A, int *ipiv,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pddesc2ge(plasma_desc_t A,
                      double *pA, int lda,
                      plasma_sequence_t *sequence,
                      plasma_request_t *request);

void plasma_pddesc2pb(plasma_desc_t A,
                      double *pA, int lda,
                      plasma_sequence_t *sequence,
                      plasma_request_t *request);

void plasma_pdge2desc(double *pA, int lda,
                      plasma_desc_t A,
                      plasma_sequence_t *sequence,
                      plasma_request_t *request);

void plasma_pdgeadd(plasma_enum_t transa,
                    double alpha,  plasma_desc_t A,
                    double beta,   plasma_desc_t B,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pdgelqf(plasma_desc_t A, plasma_desc_t T,
                    plasma_workspace_t work,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pdgelqfrh(plasma_desc_t A, plasma_desc_t T,
                      plasma_workspace_t work,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pdgemm(plasma_enum_t transa, plasma_enum_t transb,
                   double alpha, plasma_desc_t A,
                                             plasma_desc_t B,
                   double beta,  plasma_desc_t C,
                   plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pdgeqrf(plasma_desc_t A, plasma_desc_t T,
                    plasma_workspace_t work,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pdgeqrfrh(plasma_desc_t A, plasma_desc_t T,
                      plasma_workspace_t work,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pdgetri_aux(plasma_desc_t A, plasma_desc_t W,
                        plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pdgetrf(plasma_desc_t A, int *ipiv,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pdsymm(plasma_enum_t side, plasma_enum_t uplo,
                   double alpha, plasma_desc_t A,
                                             plasma_desc_t B,
                   double beta,  plasma_desc_t C,
                   plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pdsyr2k(plasma_enum_t uplo, plasma_enum_t trans,
                    double alpha, plasma_desc_t A,
                                              plasma_desc_t B,
                    double beta,              plasma_desc_t C,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pdsyrk(plasma_enum_t uplo, plasma_enum_t trans,
                   double alpha, plasma_desc_t A,
                   double beta,  plasma_desc_t C,
                   plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pdlacpy(plasma_enum_t uplo, plasma_desc_t A, plasma_desc_t B,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pdlange(plasma_enum_t norm,
                    plasma_desc_t A, double *work, double *value,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pdlansy(plasma_enum_t norm, plasma_enum_t uplo,
                    plasma_desc_t A, double *work, double *value,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pdlansy(plasma_enum_t norm, plasma_enum_t uplo,
                    plasma_desc_t A, double *work, double *value,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pdlantr(plasma_enum_t norm, plasma_enum_t uplo, plasma_enum_t diag,
                    plasma_desc_t A, double *work, double *value,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pdlascl(plasma_enum_t uplo,
                    double cfrom, double cto,
                    plasma_desc_t A,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pdlaset(plasma_enum_t uplo,
                    double alpha, double beta,
                    plasma_desc_t A,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pdlaswp(plasma_enum_t colrow,
                    plasma_desc_t A, int *ipiv, int incx,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pdlauum(plasma_enum_t uplo, plasma_desc_t A,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pdpb2desc(double *pA, int lda,
                      plasma_desc_t A,
                      plasma_sequence_t *sequence,
                      plasma_request_t *request);

void plasma_pdpbtrf(plasma_enum_t uplo, plasma_desc_t A,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pdpotrf(plasma_enum_t uplo, plasma_desc_t A,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pdsymm(plasma_enum_t side, plasma_enum_t uplo,
                   double alpha, plasma_desc_t A,
                                             plasma_desc_t B,
                   double beta,  plasma_desc_t C,
                   plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pdsyr2k(plasma_enum_t uplo, plasma_enum_t trans,
                    double alpha, plasma_desc_t A,
                                              plasma_desc_t B,
                    double beta,  plasma_desc_t C,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pdsyrk(plasma_enum_t uplo, plasma_enum_t trans,
                   double alpha, plasma_desc_t A,
                   double beta,  plasma_desc_t C,
                   plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pdtbsm(plasma_enum_t side, plasma_enum_t uplo,
                   plasma_enum_t trans, plasma_enum_t diag,
                   double alpha, plasma_desc_t A,
                                             plasma_desc_t B,
                   const int *ipiv,
                   plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pdtradd(plasma_enum_t uplo, plasma_enum_t transa,
                    double alpha,  plasma_desc_t A,
                    double beta,   plasma_desc_t B,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pdtrmm(plasma_enum_t side, plasma_enum_t uplo,
                   plasma_enum_t trans, plasma_enum_t diag,
                   double alpha, plasma_desc_t A,
                                             plasma_desc_t B,
                   plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pdtrsm(plasma_enum_t side, plasma_enum_t uplo,
                   plasma_enum_t trans, plasma_enum_t diag,
                   double alpha, plasma_desc_t A,
                                             plasma_desc_t B,
                   plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pdtrtri(plasma_enum_t uplo, plasma_enum_t diag,
                    plasma_desc_t A,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pdorglq(plasma_desc_t A, plasma_desc_t T, plasma_desc_t Q,
                    plasma_workspace_t work,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pdorglqrh(plasma_desc_t A, plasma_desc_t T, plasma_desc_t Q,
                      plasma_workspace_t work,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pdorgqr(plasma_desc_t A, plasma_desc_t T, plasma_desc_t Q,
                    plasma_workspace_t work,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pdorgqrrh(plasma_desc_t A, plasma_desc_t T, plasma_desc_t Q,
                      plasma_workspace_t work,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pdormlq(plasma_enum_t side, plasma_enum_t trans,
                    plasma_desc_t A, plasma_desc_t T, plasma_desc_t B,
                    plasma_workspace_t work,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pdormlqrh(plasma_enum_t side, plasma_enum_t trans,
                      plasma_desc_t A, plasma_desc_t T, plasma_desc_t B,
                      plasma_workspace_t work,
                      plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pdormqr(plasma_enum_t side, plasma_enum_t trans,
                    plasma_desc_t A, plasma_desc_t T, plasma_desc_t B,
                    plasma_workspace_t work,
                    plasma_sequence_t *sequence, plasma_request_t *request);

void plasma_pdormqrrh(plasma_enum_t side, plasma_enum_t trans,
                      plasma_desc_t A, plasma_desc_t T, plasma_desc_t B,
                      plasma_workspace_t work,
                      plasma_sequence_t *sequence, plasma_request_t *request);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // ICL_PLASMA_INTERNAL_D_H
