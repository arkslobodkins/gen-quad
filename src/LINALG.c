#include "LINALG.h"

#include <stdlib.h>
#include <omp.h>
#include "../plasma-17.1/include/plasma.h"


void DGEMM_LAPACK(CMatrix A, CMatrix B, CMatrix C)
{
   char TRANS1  = 'N';
   char TRANS2  = 'N';
   double alpha = 1.0;
   double beta  = 0.0;
   int M = A.rows;
   int K = A.cols;
   int N = B.cols;
   dgemm_(&TRANS1, &TRANS2, &M, &N, &K, &alpha, A.id, &M, B.id, &K, &beta, C.id, &M);
}

int DGEQR2_LAPACK(CMatrix A, Vector TAU)
{
   assert(TAU.len = MIN(A.rows, A.cols));

   int INFO;
   int LDA = A.rows;
   double *WORK = (double *)malloc(SIZE_DOUBLE(2*A.cols));
   dgeqr2_(&A.rows, &A.cols, A.id, &LDA, TAU.id, WORK, &INFO);
   free(WORK);

   return INFO;
}

int DGEQRF_LAPACK(CMatrix A, Vector TAU)
{
   assert(TAU.len = MIN(A.rows, A.cols));

   int INFO;
   int LDA = A.rows;
   int LWORK = 2*A.cols;
   double *WORK = (double *)malloc(SIZE_DOUBLE(LWORK));
   dgeqrf_(&A.rows, &A.cols, A.id, &LDA, TAU.id, WORK, &LWORK, &INFO);
   free(WORK);

   return INFO;
}

int DORMQR_LAPACK(char SIDE, char TRANS, Vector TAU, CMatrix Q, CMatrix A)
{
   int INFO;
   int LDA = A.rows;
   int LDQ = Q.rows;
   int lworkQR = SIDE == 'L' ? 2*A.cols : 2*A.rows;
   double *workQR = (double *)malloc(SIZE_DOUBLE(lworkQR));
   dormqr_(&SIDE, &TRANS, &A.rows, &A.cols, &TAU.len, Q.id, &LDQ, TAU.id, A.id, &LDA, workQR, &lworkQR, &INFO);
   free(workQR);

   return INFO;
}

int DORGQR_LAPACK(CMatrix Q, Vector TAU)
{
   int INFO;
   int LDA = Q.rows;
   int LWORK = 2*Q.cols;
   double *WORK = (double *)malloc(LWORK*sizeof(LWORK));

   dorgqr_(&Q.rows, &Q.cols, &TAU.len, Q.id, &LDA,
          TAU.id, WORK, &LWORK, &INFO);

   free(WORK);
   return INFO;
}

int DGELS_LAPACK(CMatrix A, Vector RHS_TO_X)
{
   assert(RHS_TO_X.len == MAX(A.rows, A.cols));

   int INFO;
   char TRANS = 'N';
   int NRHS = 1;
   int LDA = A.rows;
   int LEAD_DIM = MAX(A.rows, A.cols);
   int LWORK = 5*A.cols;
   Vector WORK = Vector_init(LWORK);

   dgels_(&TRANS, &A.rows, &A.cols, &NRHS, A.id, &LDA,
          RHS_TO_X.id, &LEAD_DIM, WORK.id, &WORK.len, &INFO);
   Vector_free(WORK);

   return INFO;
}

#ifdef _OPENMP
int DGELS_PLASMA(CMatrix A, Vector RHS_TO_X)
{
   assert(RHS_TO_X.len == MAX(A.rows, A.cols));

   int NRHS = 1;
   int LDA = A.rows;
   int LEAD_DIM = MAX(A.rows, A.cols);

   plasma_init();
   plasma_desc_t T;
   int INFO = plasma_dgels(PlasmaNoTrans,
                           A.rows, A.cols, NRHS,
                           A.id, LDA, &T,
                           RHS_TO_X.id, LEAD_DIM);
   plasma_desc_destroy(&T);
   plasma_finalize();

   return INFO;
}

int DGEMM_PLASMA(CMatrix A, CMatrix B, CMatrix C)
{
   char TRANS1  = PlasmaNoTrans;
   char TRANS2  = PlasmaNoTrans;
   double alpha = 1.0;
   double beta  = 0.0;
   int M = A.rows;
   int K = A.cols;
   int N = B.cols;
   int INFO = plasma_dgemm(TRANS1, TRANS2, M, N, K, alpha, A.id, M,
                           B.id, K, beta, C.id, M);
   return INFO;
}
#endif

// Simple transpose
void Transpose(int M, int N, const double *A, double *B)
{
   for(int i = 0; i < M; ++i)
      for(int j = 0; j < N; ++j)
         B[i+j*M] = A[j+i*N];
}
