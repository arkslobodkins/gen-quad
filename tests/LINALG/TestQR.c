#include <stdio.h>
#include <stdlib.h>
#include <mkl.h>

#include "../../src/Vector.h"
#include "../../src/Matrix.h"
#include "../../src/LINALG.h"

CMatrix read_matrix_A()
{
   int rows, cols;
   FILE *fp = fopen("A1.txt", "r");
   fscanf(fp, "%i", &rows);
   fscanf(fp, "%i", &cols);
   CMatrix A = CMatrix_init(rows, cols);

   for(int j = 0; j < A.cols; ++j)
      for(int i = 0; i < A.rows; ++i)
         fscanf(fp, "%lf", &A.cid[j][i]);

   fclose(fp);
   return A;
}

CMatrix read_matrix_Q()
{
   FILE *fp = fopen("Q1.txt", "r");
   int rows, cols;
   fscanf(fp, "%i", &rows);
   fscanf(fp, "%i", &cols);
   CMatrix Q = CMatrix_init(rows, cols);

   for(int j = 0; j < Q.cols; ++j)
      for(int i = 0; i < Q.rows; ++i)
         fscanf(fp, "%lf", &Q.cid[j][i]);

   fclose(fp);
   return Q;
}

int main(int argc, char *argv[])
{
   mkl_set_threading_layer(MKL_THREADING_SEQUENTIAL);

   CMatrix A = read_matrix_A(A);
   CMatrix Q_MATLAB = read_matrix_Q(Q_MATLAB);
   CMatrix Q_LINALG = CMatrix_init(Q_MATLAB.rows, Q_MATLAB.cols);
   Vector TAU = Vector_init(MIN(A.rows, A.cols));

   DGEQRF_LAPACK(A, TAU);
   CMatrix_Identity(Q_LINALG);
   DORMQR_LAPACK('L', 'N', TAU, A, Q_LINALG);

   printf("Maximum difference = %.4e\n", CMatrixMaxDifference(Q_MATLAB, Q_LINALG));

   CMatrix_free(A);
   CMatrix_free(Q_LINALG);
   CMatrix_free(Q_MATLAB);
   Vector_free(TAU);

   return EXIT_SUCCESS;
}
