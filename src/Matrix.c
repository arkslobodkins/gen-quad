/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */

#include "Matrix.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

Matrix Matrix_init(int nRows, int nCols)
{
   Matrix M = {0};

   M.rows = nRows;
   M.cols = nCols;
   M.id = (double **)calloc(nRows, sizeof(double *));

   for(int i = 0; i < nRows; ++i)
      M.id[i] = (double *)calloc(nCols, sizeof(double *));

   return M;
}


void Matrix_realloc(int nRows, int nCols, Matrix *M)
{
   M->rows = nRows;
   M->cols = nCols;
   M->id = (double **)realloc(M->id, nRows* sizeof(double *));

   for(int i = 0; i < nRows; ++i)
      M->id[i] = (double *)realloc(M->id[i], nCols * sizeof(double *));
}


void Matrix_free(Matrix M)
{
   if(M.id != NULL)
   {
      for(int i = 0; i < M.rows; ++i) free(M.id[i]);
      free(M.id); M.id = NULL;
   }
}


void MatVec(Matrix M, Vector V, Vector O)
{
   assert(M.cols == V.len);
   assert(O.len >= M.rows);

   int i, j;

   for(i = 0; i < M.rows; ++i)
      for(j = 0; j < M.cols; ++j)
         O.id[i] = 0.0;

   for(i = 0; i < M.rows; ++i)
      for(j = 0; j < M.cols; ++j)
         O.id[i] += M.id[i][j] * V.id[j];

}


void PrintMatrix(const Matrix M)
{
   for(int i = 0; i < M.rows; ++i)
   {
      for(int j = 0; j < M.cols; ++j)
         printf("M[%i][%i] = %f \n", i, j, M.id[i][j]);

      printf("\n");
   }
}
