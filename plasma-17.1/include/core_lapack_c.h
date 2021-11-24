/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from include/core_lapack_z.h, normal z -> c, Mon Nov 22 19:09:34 2021
 *
 **/

#ifndef ICL_CORE_LAPACK_C_H
#define ICL_CORE_LAPACK_C_H

#ifdef __cplusplus
extern "C" {
#endif

// LAPACK_GLOBAL is Fortran name mangling macro from LAPACKE

// LAPACKE_clantr broken (returns 0) in LAPACKE < 3.6.1
#ifndef LAPACK_clantr
#define LAPACK_clantr LAPACK_GLOBAL(clantr, CLANTR)
float LAPACK_clantr(const char *norm, const char *uplo, const char *diag,
                     const int *m, const int *n,
                     const plasma_complex32_t *A, const int *lda,
                     float *work);
#endif

// LAPACKE_clascl not available in LAPACKE < 3.6.0
#ifndef LAPACK_clascl
#define LAPACK_clascl LAPACK_GLOBAL(clascl, CLASCL)
void LAPACK_clascl(const char *type, const int *kl, const int *ku,
                   const float *cfrom, const float *cto,
                   const int *m, const int *n,
                   plasma_complex32_t *A, const int *lda,
                   int *info);
#endif

// LAPACKE_classq not available yet
#ifndef LAPACK_classq
#define LAPACK_classq LAPACK_GLOBAL(classq, CLASSQ)
void LAPACK_classq(const int *n, const plasma_complex32_t *x, const int *incx,
                   float *scale, float *sumsq);
#endif

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // ICL_CORE_LAPACK_C_H
