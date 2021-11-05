#include <stdlib.h>

#include "LINALG.h"
#include "GENERAL_QUADRATURE.h"
#include "Matrix.h"


void MAT_MUL(int M, int K, int N, double *MAT1, double *MAT2, double *MAT3)
{
   char TRANS1  = 'N';
   char TRANS2  = 'N';
   double alpha = 1.0;
   double beta  = 0.0;
   dgemm_(&TRANS1, &TRANS2, &M, &N, &K, &alpha, MAT1, &M, MAT2, &K, &beta, MAT3, &M);

}


int DORMQR_M(char SIDE, char TRANS, double *REFL_Q, CMatrix Q, CMatrix A)
{
   int INFO;
   int LDJ        = Q.rows;
   int N_REFL     = MIN(Q.rows, Q.cols);
   int lworkQR    = 2*Q.cols;
   double *workQR = (double *) malloc( SIZE_DOUBLE(lworkQR));
   dormqr_(&SIDE, &TRANS, &A.rows, &A.cols, &N_REFL, Q.id, &LDJ, REFL_Q, A.id, &LDJ, workQR, &lworkQR, &INFO);

   free(workQR);
   return INFO;
}


int DORMQR_V(char SIDE, char TRANS, double *REFL_Q, CMatrix Q, Vector A)
{
   int INFO;
   int LDJ        = Q.rows;
   int ONE_COLUMN = 1;
   int N_REFL     = MIN(Q.rows, Q.cols);
   int lworkQR    = 2*Q.cols;
   double *workQR = (double *) malloc( SIZE_DOUBLE(lworkQR));
   dormqr_(&SIDE, &TRANS, &A.len, &ONE_COLUMN, &N_REFL, Q.id, &LDJ, REFL_Q, A.id, &LDJ, workQR, &lworkQR, &INFO);

   free(workQR);
   return INFO;
}


void Transpose(int M, int N, const double *A, double *B)
{
   for(int i = 0; i < M; ++i)
      for(int j = 0; j < N; ++j)
         B[i+j*M] = A[j+i*N];
}
