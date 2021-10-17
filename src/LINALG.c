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


void C_TO_FORTRAN(int M, int N, const double *A, double *B)
{
   for(int i = 0; i < M; ++i)
      for(int j = 0; j < N; ++j)
         B[ ij2(j,i,M) ] = A[ ij2(i,j,N) ];

}


void FORTRAN_TO_C(int M, int N, const double *A, double *B)
{
   for(int i = 0; i < M; ++i)
      for(int j = 0; j < N; ++j)
         B[ ij2(i,j,N) ] = A[ ij2(j,i,M) ];

}


void MATRIX_TO_FORTRAN(int M, int N, const Matrix A, double *B)
{
   for(int i = 0; i < M; ++i)
      for(int j = 0; j < N; ++j)
         B[ ij2(j,i,M) ] = A.id[i][j];

}


void FORTRAN_TO_MATRIX(int M, int N, const double *A, Matrix B)
{
   for(int i = 0; i < M; ++i)
      for(int j = 0; j < N; ++j)
         B.id[i][j] = A[ ij2(j,i,M) ];

}
