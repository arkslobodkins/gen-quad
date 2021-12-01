/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */

#include "Matrix.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <omp.h>

RMatrix RMatrix_init(int nRows, int nCols)
{
   RMatrix M = {0};

   M.rows = nRows;
   M.cols = nCols;
   M.len  = nRows*nCols;
   M.id   = (double *)calloc(nRows*nCols, sizeof(double));
   M.rid  = (double **)malloc(nRows*sizeof(double*));

   for(int i = 0; i < nRows; ++i)
      M.rid[i] = &M.id[nCols*i];

   return M;
}

void RMatrix_realloc(int nRows, int nCols, RMatrix *M)
{
   M->rows = nRows;
   M->cols = nCols;
   M->len  = nRows*nCols;

   M->rid = (double **)realloc(M->rid, nRows*sizeof(double*));
   M->id  = (double *)realloc(M->id, nRows*nCols*sizeof(double));
   memset(M->id, 0, M->len*sizeof(double));

   for(int i = 0; i < nRows; ++i)
      M->rid[i] = &M->id[nCols*i];

}

void RMatrix_free(RMatrix M)
{
   if(M.id != NULL)  { free(M.id); M.id = NULL; }
   if(M.rid != NULL) { free(M.rid); M.rid = NULL; }
}

void RMatVec(RMatrix M, Vector x, Vector y)
{
   assert(M.cols == x.len);
   assert(y.len == M.rows);

   memset(y.id, 0, y.len*sizeof(double));
   for(int i = 0; i < M.rows; ++i)
      for(int j = 0; j < M.cols; ++j)
         y.id[i] += R_ELEM_ID(M, i, j) * x.id[j];
}

void RMatrixPrint(RMatrix M)
{
   for(int i = 0; i < M.rows; ++i)
   {
      for(int j = 0; j < M.cols; ++j)
         printf("M[%i][%i] = %.10e \n", i, j, R_ELEM_ID(M, i, j));

      printf("\n");
   }
}



CMatrix CMatrix_init(int nRows, int nCols)
{
   CMatrix M = {0};

   M.rows = nRows;
   M.cols = nCols;
   M.len  = nRows*nCols;
   M.id   = (double *)calloc(nRows*nCols, sizeof(double));
   M.cid  = (double **)malloc(nCols*sizeof(double*));

   for(int j = 0; j < nCols; ++j)
      M.cid[j] = &M.id[nRows*j];
   return M;
}

void CMatrix_realloc(int nRows, int nCols, CMatrix *M)
{
   M->rows = nRows;
   M->cols = nCols;
   M->len  = nRows*nCols;

   M->id  = (double *)realloc(M->id, nRows*nCols*sizeof(double));
   M->cid = (double **)realloc(M->cid, nCols*sizeof(double*));
   memset(M->id, 0, M->len*sizeof(double));

   for(int j = 0; j < nCols; ++j)
      M->cid[j] = &M->id[nRows*j];
}

void CMatrix_free(CMatrix M)
{
   if(M.id != NULL)  { free(M.id); M.id = NULL; }
   if(M.cid != NULL) { free(M.cid); M.cid = NULL; }
}

void CMatrix_Assign(CMatrix A, CMatrix B)
{
   assert(A.rows == B.rows);
   assert(A.cols == B.cols);
   memcpy(B.id, A.id, A.len*sizeof(double));
}

CMatrix CMatrix_Transpose(CMatrix A)
{
   CMatrix A_TR = CMatrix_init(A.cols, A.rows);
   int Acols = A.cols;
   int Arows = A.rows;

   #ifdef _OPENMP
   #pragma omp parallel for default(shared) schedule(static) num_threads(omp_get_max_threads())
   #endif
   for(int j = 0; j < Acols; ++j)
      for(int i = 0; i < Arows; ++i)
         A_TR.cid[i][j] = A.cid[j][i];
   CMatrix_free(A);

   return A_TR;
}

void CMatVec(CMatrix M, Vector x, Vector y)
{
   assert(M.cols == x.len);
   assert(y.len == M.rows);

   memset(y.id, 0, y.len*sizeof(double));
   for(int j = 0; j < M.cols; ++j)
      for(int i = 0; i < M.rows; ++i)
         y.id[i] += C_ELEM_ID(M, i, j) * x.id[j];
}

void CMatrixPrint(CMatrix M)
{
   for(int i = 0; i < M.rows; ++i)
   {
      for(int j = 0; j < M.cols; ++j)
         printf("M[%i][%i] = %.10e \n", i, j, C_ELEM_ID(M, i, j));

      printf("\n");
   }
}
