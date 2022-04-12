/* Arkadijs Slobodkins
 * SMU Mathematics
 * Copyright 2022
 */


#include "LINALG.h"
#include "Print.h"

#include <stdlib.h>
#include <omp.h>
#include <plasma.h>
#include <mkl_blas.h>
#include <mkl_lapacke.h>

void DGEMM_LAPACK(CMatrix A, CMatrix B, CMatrix C)
{
   assert(C.rows == A.rows);
   assert(C.cols == B.cols);
   assert(A.cols == B.rows);

   char TRANS1  = 'N';
   char TRANS2  = 'N';
   double alpha = 1.0;
   double beta  = 0.0;
   int M = A.rows;
   int K = A.cols;
   int N = B.cols;

   dgemm( &TRANS1, &TRANS2, &M, &N, &K, &alpha, A.id, &M, B.id, &K, &beta, C.id, &M );
}

int DGEQR2_LAPACK(CMatrix A, Vector TAU)
{
   assert(TAU.len == MIN(A.rows, A.cols));

   int LDA = A.rows;
   int INFO = LAPACKE_dgeqr2( LAPACK_COL_MAJOR, A.rows, A.cols, A.id, LDA, TAU.id );

   if(INFO != 0) return LAPACK_ERR;
   return GQ_SUCCESS;
}

int DGEQRF_LAPACK(CMatrix A, Vector TAU)
{
   assert(TAU.len == MIN(A.rows, A.cols));

   int LDA = A.rows;
   int INFO = LAPACKE_dgeqrf( LAPACK_COL_MAJOR, A.rows, A.cols, A.id, LDA, TAU.id );

   if(INFO != 0) return LAPACK_ERR;
   return GQ_SUCCESS;
}

int DORMQR_LAPACK(char SIDE, char TRANS, Vector TAU, CMatrix Q, CMatrix A)
{
   assert(TAU.len == MIN(Q.rows, Q.cols));

   int LDA;
   if(SIDE == 'L')      LDA = A.rows;
   else if(SIDE == 'R') LDA = A.cols;
   else                 return INV_INPUT;
   int LDQ = Q.rows;

   int INFO = LAPACKE_dormqr( LAPACK_COL_MAJOR, SIDE, TRANS,
                              A.rows, A.cols, TAU.len,
                              Q.id, LDQ, TAU.id, A.id, LDA );

   if(INFO != 0) return LAPACK_ERR;
   return GQ_SUCCESS;
}

int DORGQR_LAPACK(CMatrix Q, Vector TAU)
{
   assert(TAU.len == MIN(Q.rows, Q.cols));

   int LDA = Q.rows;
   int INFO = LAPACKE_dorgqr( LAPACK_COL_MAJOR, Q.rows, Q.cols,
                              TAU.len, Q.id, LDA, TAU.id );
   if(INFO != 0) return LAPACK_ERR;
   return GQ_SUCCESS;
}

int DGESVD_LAPACK(CMatrix A, Vector SINGV)
{
   assert(SINGV.len = MIN(A.rows, A.cols));

   int M = A.rows;
   int N = A.cols;
   int LDA = A.rows;

   char JOBU  = 'N'; // U is not computed
   char JOBVT = 'N'; // V is not computed
   double *U  = NULL;
   double *VT = NULL;
   int LDU  = 1;
   int LDVT = 1;
   Vector SUPERB = Vector_init(MIN(M, N)-1);

   int INFO = LAPACKE_dgesvd( LAPACK_COL_MAJOR, JOBU, JOBVT,
                               M, N, A.id, LDA, SINGV.id, U, LDU,
                               VT, LDVT, SUPERB.id );
   Vector_free(SUPERB);

   if(INFO != 0) return LAPACK_ERR;
   return GQ_SUCCESS;
}

int DGELS_LAPACK(CMatrix A, Vector b, Vector x)
{
   assert(A.rows == b.len);
   assert(A.cols == x.len);

   char TRANS = 'N';
   int rows = A.rows;
   int cols = A.cols;
   int NRHS = 1;
   int LDA = A.rows;
   int LEAD_DIM = MAX(A.rows, A.cols);

   Vector RHS_TO_X = Vector_init(MAX(rows, cols));
   for(int i = 0; i < rows; ++i)
      RHS_TO_X.id[i] = b.id[i];

   int INFO = LAPACKE_dgels( LAPACK_COL_MAJOR, TRANS, A.rows,
                             A.cols, NRHS, A.id, LDA, RHS_TO_X.id, LEAD_DIM );

   for(int i = 0; i < cols; ++i)
      x.id[i] = RHS_TO_X.id[i];
   Vector_free(RHS_TO_X);

   if(INFO > 0)
      PRINT_WARN("LAPACK encountered 0 pivot", __LINE__, __FILE__);

   if(INFO != 0)
      return LAPACK_ERR;
   return GQ_SUCCESS;
}

double MIN_SINGV_LAPACK(CMatrix A)
{
   Vector singv = Vector_init(MIN(A.rows, A.cols));

   int INFO = DGESVD_LAPACK(A, singv);
   if(INFO != GQ_SUCCESS)
   {
      PRINT_ERR(STR_LAPACK_ERR, __LINE__, __FILE__);
      return -1.;
   }

   VMin min_singv = VectorMin(singv);

   Vector_free(singv);
   return min_singv.min_value;
}

#ifdef _OPENMP

int DGEMM_PLASMA(CMatrix A, CMatrix B, CMatrix C)
{
   assert(C.rows == A.rows);
   assert(C.cols == B.cols);
   assert(A.cols == B.rows);

   char TRANS1 = PlasmaNoTrans;
   char TRANS2 = PlasmaNoTrans;
   int M = A.rows;
   int K = A.cols;
   int N = B.cols;
   double alpha = 1.0;
   double beta  = 0.0;

   plasma_init();
   int INFO = plasma_dgemm( TRANS1, TRANS2, M, N, K, alpha, A.id, M,
                            B.id, K, beta, C.id, M );
   plasma_finalize();

   if(INFO != 0) return PLASMA_ERR;
   return GQ_SUCCESS;
}

int DGELS_PLASMA(CMatrix A, Vector b, Vector x)
{
   assert(A.rows == b.len);
   assert(A.cols == x.len);

   int rows = A.rows;
   int cols = A.cols;
   int NRHS = 1;
   int LDA = A.rows;
   int LEAD_DIM = MAX(A.rows, A.cols);

   Vector RHS_TO_X = Vector_init(MAX(rows, cols));
   for(int i = 0; i < rows; ++i)
      RHS_TO_X.id[i] = b.id[i];

   plasma_init();
   plasma_desc_t T;
   int INFO = plasma_dgels( PlasmaNoTrans,
                            A.rows, A.cols, NRHS,
                            A.id, LDA, &T,
                            RHS_TO_X.id, LEAD_DIM );
   plasma_desc_destroy(&T);
   plasma_finalize();

   for(int i = 0; i < cols; ++i)
      x.id[i] = RHS_TO_X.id[i];
   Vector_free(RHS_TO_X);

   if(INFO != 0) return PLASMA_ERR;
   return GQ_SUCCESS;
}

#endif

