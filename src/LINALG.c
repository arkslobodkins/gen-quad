/* Arkadijs Slobodkins
 * SMU Mathematics
 * Copyright 2022
 */


#include "Print.h"

#include <stdlib.h>
#include <omp.h>
#include <plasma.h>
#include <mkl_blas.h>
#include <mkl_lapacke.h>

int DGEMM_LAPACK(CMatrix A, CMatrix B, CMatrix C)
{
   if(A.rows != C.rows) return INV_INPUT;
   if(A.cols != B.rows) return INV_INPUT;
   if(B.cols != C.cols) return INV_INPUT;

   char TRANS1  = 'N';
   char TRANS2  = 'N';
   double alpha = 1.0;
   double beta  = 0.0;
   int M = A.rows;
   int K = A.cols;
   int N = B.cols;

   dgemm( &TRANS1, &TRANS2, &M, &N, &K, &alpha, A.id, &M, B.id, &K, &beta, C.id, &M );
   return GQ_SUCCESS;
}

int DGEQR2_LAPACK(CMatrix A, Vector TAU)
{
   assert(TAU.len == MIN(A.rows, A.cols));

   int LDA = A.rows;
   int INFO = LAPACKE_dgeqr2( LAPACK_COL_MAJOR, A.rows, A.cols, A.id, LDA, TAU.id );

   if(INFO != 0) return LAPACK_ERR;
   else          return GQ_SUCCESS;
}

int DGEQRF_LAPACK(CMatrix A, Vector TAU)
{
   assert(TAU.len == MIN(A.rows, A.cols));

   int LDA = A.rows;
   int INFO = LAPACKE_dgeqrf( LAPACK_COL_MAJOR, A.rows, A.cols, A.id, LDA, TAU.id );

   if(INFO != 0) return LAPACK_ERR;
   else          return GQ_SUCCESS;
}

int DORMQR_LAPACK(char SIDE, char TRANS, Vector TAU, CMatrix Q, CMatrix A)
{
   assert(TAU.len == MIN(Q.rows, Q.cols));

   int LDA;
   if(SIDE == 'L') LDA = A.rows;
   else if(SIDE == 'R') LDA = A.cols;
   else            return INV_INPUT;

   int LDQ = Q.rows;
   int INFO = LAPACKE_dormqr( LAPACK_COL_MAJOR, SIDE, TRANS,
                              A.rows, A.cols, TAU.len,
                              Q.id, LDQ, TAU.id, A.id, LDA );

   if(INFO != 0) return LAPACK_ERR;
   else          return GQ_SUCCESS;
}

int DORGQR_LAPACK(CMatrix Q, Vector TAU)
{
   int LDA = Q.rows;
   int INFO = LAPACKE_dorgqr( LAPACK_COL_MAJOR, Q.rows, Q.cols,
                              TAU.len, Q.id, LDA, TAU.id );
   if(INFO != 0) return LAPACK_ERR;
   else          return GQ_SUCCESS;
}

int DGESVD_LAPACK(CMatrix A, CMatrix VT)
{
   assert(VT.rows == MIN(A.rows, A.cols));
   assert(VT.cols == A.cols);
   char JOBU  = 'N'; // U is not computed
   char JOBVT = 'S'; //  return MIN(M, N) rows of V^T, i.e. right singular vectors

   int M = A.rows;
   int N = A.cols;
   int LDA = A.rows;
   Vector SINGV = Vector_init(MIN(M, N));
   if(N > M)
      PRINT_WARN("svd received underdetermined matrix", __LINE__, __FILE__);

   double *U = NULL;
   int LDU = 1;
   int LDVT = MIN(M, N);    // expected to be N in general
   int LWORK = 8*MIN(M, N); // assumes M is larger than N, 5*MIN is the minimum work required
   Vector WORK = Vector_init(LWORK);
   int INFO = LAPACKE_dgesvd( LAPACK_COL_MAJOR, JOBU, JOBVT,
                               M, N, A.id, LDA, SINGV.id, U, LDU,
                               VT.id, LDVT, WORK.id );
   Vector_free(SINGV);
   Vector_free(WORK);

   if(INFO != 0) return LAPACK_ERR;
   else          return GQ_SUCCESS;
}

int DGELS_LAPACK(CMatrix A, Vector b, Vector x)
{
   assert(A.rows == b.len);
   assert(A.cols == x.len);
   int rows = A.rows;
   int cols = A.cols;

   Vector RHS_TO_X = Vector_init(MAX(rows, cols));
   for(int i = 0; i < rows; ++i) RHS_TO_X.id[i] = b.id[i];

   char TRANS = 'N';
   int NRHS = 1;
   int LDA = A.rows;
   int LEAD_DIM = MAX(A.rows, A.cols);
   int INFO = LAPACKE_dgels( LAPACK_COL_MAJOR, TRANS, A.rows,
                             A.cols, NRHS, A.id, LDA, RHS_TO_X.id, LEAD_DIM );

   for(int i = 0; i < cols; ++i) x.id[i] = RHS_TO_X.id[i];
   Vector_free(RHS_TO_X);

   if(INFO > 0)
      PRINT_WARN("LAPACK encountered 0 pivot", __LINE__, __FILE__);

   if(INFO != 0) return LAPACK_ERR;
   else          return GQ_SUCCESS;

}

//int DGELS_LAPACK(CMatrix A, Vector RHS_TO_X)
//{
//   assert(RHS_TO_X.len == MAX(A.rows, A.cols));
//
//   char TRANS = 'N';
//   int NRHS = 1;
//   int LDA = A.rows;
//   int LEAD_DIM = MAX(A.rows, A.cols);
//   int INFO = LAPACKE_dgels( LAPACK_COL_MAJOR, TRANS, A.rows,
//                             A.cols, NRHS, A.id, LDA, RHS_TO_X.id, LEAD_DIM );
//   if(INFO > 0)
//      PRINT_WARN("LAPACK encountered 0 pivot", __LINE__, __FILE__);
//
//   if(INFO != 0) return LAPACK_ERR;
//   else          return GQ_SUCCESS;
//
//}

#ifdef _OPENMP
int DGELS_PLASMA(CMatrix A, Vector b, Vector x)
{
   assert(A.rows == b.len);
   assert(A.cols == x.len);
   int rows = A.rows;
   int cols = A.cols;

   Vector RHS_TO_X = Vector_init(MAX(rows, cols));
   for(int i = 0; i < rows; ++i) RHS_TO_X.id[i] = b.id[i];

   int NRHS = 1;
   int LDA = A.rows;
   int LEAD_DIM = MAX(A.rows, A.cols);

   plasma_init();
   plasma_desc_t T;
   int INFO = plasma_dgels( PlasmaNoTrans,
                            A.rows, A.cols, NRHS,
                            A.id, LDA, &T,
                            RHS_TO_X.id, LEAD_DIM );
   plasma_desc_destroy(&T);
   plasma_finalize();

   for(int i = 0; i < cols; ++i) x.id[i] = RHS_TO_X.id[i];
   Vector_free(RHS_TO_X);

   if(INFO != 0) return PLASMA_ERR;
   else          return GQ_SUCCESS;
}

int DGEMM_PLASMA(CMatrix A, CMatrix B, CMatrix C)
{
   if(A.rows != C.rows) return INV_INPUT;
   if(A.cols != B.rows) return INV_INPUT;
   if(B.cols != C.cols) return INV_INPUT;

   plasma_init();
   char TRANS1  = PlasmaNoTrans;
   char TRANS2  = PlasmaNoTrans;
   double alpha = 1.0;
   double beta  = 0.0;
   int M = A.rows;
   int K = A.cols;
   int N = B.cols;
   int INFO = plasma_dgemm( TRANS1, TRANS2, M, N, K, alpha, A.id, M,
                            B.id, K, beta, C.id, M );
   plasma_finalize();
   if(INFO != 0) return PLASMA_ERR;
   else          return GQ_SUCCESS;
}

#endif

void Transpose(int M, int N, const double *A, double *B)
{
   for(int i = 0; i < M; ++i)
      for(int j = 0; j < N; ++j)
         B[i+j*M] = A[j+i*N];
}
