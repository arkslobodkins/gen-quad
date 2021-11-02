#include "LINALG.h"
#include "GENERAL_QUADRATURE.h"
#include "Matrix.h"


void MAT_MUL(int M, int K, int N, double *MAT1, double *MAT2, double *MAT3)
{
   char TRANS1 = 'N';
   char TRANS2 = 'N';
   double alpha = 1.0;
   double beta = 0.0;
   dgemm_(&TRANS1, &TRANS2, &M, &N, &K, &alpha, MAT1, &M, MAT2, &K, &beta, MAT3, &M);

}


void Transpose(int M, int N, const double *A, double *B)
{
   for(int i = 0; i < M; ++i)
      for(int j = 0; j < N; ++j)
         B[i+j*M] = A[j+i*N];
}
