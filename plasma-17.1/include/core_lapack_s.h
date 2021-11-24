/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from include/core_lapack_z.h, normal z -> s, Mon Nov 22 19:09:33 2021
 *
 **/

#ifndef ICL_CORE_LAPACK_S_H
#define ICL_CORE_LAPACK_S_H

#ifdef __cplusplus
extern "C" {
#endif

// LAPACK_GLOBAL is Fortran name mangling macro from LAPACKE

// LAPACKE_slantr broken (returns 0) in LAPACKE < 3.6.1
#ifndef LAPACK_slantr
#define LAPACK_slantr LAPACK_GLOBAL(slantr, SLANTR)
float LAPACK_slantr(const char *norm, const char *uplo, const char *diag,
                     const int *m, const int *n,
                     const float *A, const int *lda,
                     float *work);
#endif

// LAPACKE_slascl not available in LAPACKE < 3.6.0
#ifndef LAPACK_slascl
#define LAPACK_slascl LAPACK_GLOBAL(slascl, SLASCL)
void LAPACK_slascl(const char *type, const int *kl, const int *ku,
                   const float *cfrom, const float *cto,
                   const int *m, const int *n,
                   float *A, const int *lda,
                   int *info);
#endif

// LAPACKE_slassq not available yet
#ifndef LAPACK_slassq
#define LAPACK_slassq LAPACK_GLOBAL(slassq, SLASSQ)
void LAPACK_slassq(const int *n, const float *x, const int *incx,
                   float *scale, float *sumsq);
#endif

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // ICL_CORE_LAPACK_S_H
