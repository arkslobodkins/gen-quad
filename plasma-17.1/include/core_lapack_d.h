/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from include/core_lapack_z.h, normal z -> d, Mon Nov 22 19:09:33 2021
 *
 **/

#ifndef ICL_CORE_LAPACK_D_H
#define ICL_CORE_LAPACK_D_H

#ifdef __cplusplus
extern "C" {
#endif

// LAPACK_GLOBAL is Fortran name mangling macro from LAPACKE

// LAPACKE_dlantr broken (returns 0) in LAPACKE < 3.6.1
#ifndef LAPACK_dlantr
#define LAPACK_dlantr LAPACK_GLOBAL(dlantr, DLANTR)
double LAPACK_dlantr(const char *norm, const char *uplo, const char *diag,
                     const int *m, const int *n,
                     const double *A, const int *lda,
                     double *work);
#endif

// LAPACKE_dlascl not available in LAPACKE < 3.6.0
#ifndef LAPACK_dlascl
#define LAPACK_dlascl LAPACK_GLOBAL(dlascl, DLASCL)
void LAPACK_dlascl(const char *type, const int *kl, const int *ku,
                   const double *cfrom, const double *cto,
                   const int *m, const int *n,
                   double *A, const int *lda,
                   int *info);
#endif

// LAPACKE_dlassq not available yet
#ifndef LAPACK_dlassq
#define LAPACK_dlassq LAPACK_GLOBAL(dlassq, DLASSQ)
void LAPACK_dlassq(const int *n, const double *x, const int *incx,
                   double *scale, double *sumsq);
#endif

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // ICL_CORE_LAPACK_D_H
