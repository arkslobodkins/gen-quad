/* Arkadijs Slobodkins
 * SMU Mathematics
 * Copyright 2022
 */


#include "Matrix.h"
#include "Print.h"
#include "GENERAL_QUADRATURE.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mkl_blas.h>

RMatrix RMatrix_init(int nRows, int nCols)
{
   assert(nRows >= 1);
   assert(nCols >= 1);

   RMatrix M;

   M.rows = nRows;
   M.cols = nCols;
   M.len  = nRows*nCols;

   M.id = (double *)calloc(nRows*nCols, sizeof(double));
   if(M.id == NULL) PRINT_ERR(STR_ALLOC_FAIL, __LINE__, __FILE__);

   M.rid  = (double **)malloc(nRows*sizeof(double*));
   if(M.rid == NULL) PRINT_ERR(STR_ALLOC_FAIL, __LINE__, __FILE__);

   for(int i = 0; i < nRows; ++i)
      M.rid[i] = &M.id[nCols*i];

   return M;
}

void RMatrix_realloc(int nRows, int nCols, RMatrix *M)
{
   assert(nRows >= 1);
   assert(nCols >= 1);

   M->rows = nRows;
   M->cols = nCols;
   M->len  = nRows*nCols;

   M->rid = (double **)realloc(M->rid, nRows*sizeof(double*));
   if(M->rid == NULL) PRINT_ERR(STR_ALLOC_FAIL, __LINE__, __FILE__);

   M->id  = (double *)realloc(M->id, nRows*nCols*sizeof(double));
   if(M->id == NULL) PRINT_ERR(STR_ALLOC_FAIL, __LINE__, __FILE__);

   memset(M->id, 0, M->len*sizeof(double));

   for(int i = 0; i < nRows; ++i)
      M->rid[i] = &M->id[nCols*i];
}

void RMatrix_free(RMatrix M)
{
   if(M.id != NULL)  { free(M.id); }
   if(M.rid != NULL) { free(M.rid); }
}

void RMatrix_SetZero(RMatrix M)
{
   assert(M.rows > 0 && M.cols > 0);
   memset(M.id, 0, SIZE_DOUBLE(M.rows*M.cols));
}

void RMatrix_LoadToRow(int row, RMatrix M, Vector v)
{
   assert(row > -1 && row < M.rows);
   assert(v.len == M.cols);
   memcpy(M.rid[row], v.id, SIZE_DOUBLE(M.cols));
}

void RMatVec(RMatrix M, Vector x, Vector y)
{
   assert(M.rows == y.len);
   assert(M.cols == x.len);

   memset(y.id, 0, y.len*sizeof(double));
   for(int i = 0; i < M.rows; ++i)
      for(int j = 0; j < M.cols; ++j)
         y.id[i] += R_ELEM_ID(M, i, j) * x.id[j];
}

double RDotRow(int row, RMatrix M, Vector x)
{
   assert(row > -1 && row < M.rows);
   assert(M.cols == x.len);

   double dot = 0.;
   for(int j = 0; j < M.cols; ++j)
      dot += M.rid[row][j] * x.id[j];
   return dot;
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
   assert(nRows >= 1);
   assert(nCols >= 1);

   CMatrix M;

   M.rows = nRows;
   M.cols = nCols;
   M.len  = nRows*nCols;

   M.id = (double *)calloc(nRows*nCols, sizeof(double));
   if(M.id == NULL) PRINT_ERR(STR_ALLOC_FAIL, __LINE__, __FILE__);

   M.cid = (double **)malloc(nCols*sizeof(double*));
   if(M.cid == NULL) PRINT_ERR(STR_ALLOC_FAIL, __LINE__, __FILE__);

   for(int j = 0; j < nCols; ++j)
      M.cid[j] = &M.id[nRows*j];
   return M;
}

void CMatrix_realloc(int nRows, int nCols, CMatrix *M)
{
   assert(nRows >= 1);
   assert(nCols >= 1);

   M->rows = nRows;
   M->cols = nCols;
   M->len  = nRows*nCols;

   M->id = (double *)realloc(M->id, nRows*nCols*sizeof(double));
   if(M->id == NULL) PRINT_ERR(STR_ALLOC_FAIL, __LINE__, __FILE__);

   M->cid = (double **)realloc(M->cid, nCols*sizeof(double*));
   if(M->cid == NULL) PRINT_ERR(STR_ALLOC_FAIL, __LINE__, __FILE__);
   memset(M->id, 0, M->len*sizeof(double));

   for(int j = 0; j < nCols; ++j)
      M->cid[j] = &M->id[nRows*j];
}

void CMatrix_free(CMatrix M)
{
   if(M.id != NULL)  { free(M.id); }
   if(M.cid != NULL) { free(M.cid); }
}

void CMatrix_SetZero(CMatrix M)
{
   assert(M.rows > 0 && M.cols > 0);
   memset(M.id, 0, SIZE_DOUBLE(M.rows*M.cols));
}

void CMatrix_Assign(CMatrix A, CMatrix B)
{
   assert(A.rows == B.rows);
   assert(A.cols == B.cols);
   memcpy(B.id, A.id, A.len*sizeof(double));
}

void CMatrix_Identity(CMatrix M)
{
   assert(M.rows == M.cols);
   memset(M.id, 0, M.len*sizeof(double));
   for(int i = 0; i < M.rows; ++i)
      C_ELEM_ID(M, i, i) = 1.0;
}

// Load all entries from Vector x to col column of Matrix M
void CMatrix_LoadToColumn(int col, CMatrix M, Vector x)
{
   assert(col > -1 && col < M.cols);
   assert(x.len == M.rows);
   memcpy(M.cid[col], x.id, SIZE_DOUBLE(M.rows));
}

// Load entries[begin, begin+x.len-1] from Vector x to col column of Matrix M
void CMatrix_LoadToColumnRange(int begin, int col, CMatrix M, Vector x)
{
   assert(col > -1 && col < M.cols);
   assert(begin > -1);              // first entry copied from vector
   assert(x.len % M.rows == 0);     // length of x multiple of #rows
   assert(begin + M.rows <= x.len); // make sure begin is not too far

   memcpy(M.cid[col], x.id+begin, SIZE_DOUBLE(M.rows));
}

void CMatrix_GetRow(int row, CMatrix M, Vector x)
{
   assert(row > -1 && row < M.rows);
   assert(x.len == M.cols);
   for(int j = 0; j < M.cols; ++j)
      x.id[j] = C_ELEM_ID(M, row, j);
}

void CMatrix_GetColumn(int col, CMatrix M, Vector x)
{
   assert(col > -1 && col < M.cols);
   assert(x.len == M.rows);
   memcpy(x.id, M.cid[col], SIZE_DOUBLE(M.rows));
}

void CMatrix_ColumnScale(int col, double c, CMatrix M)
{
   assert(col > -1 && col < M.cols);
   int spacing = 1;
   dscal(&M.rows, &c, M.cid[col], &spacing);
}

void CMatrix_Assign_Transpose(CMatrix A, CMatrix B)
{
   assert(A.rows == B.cols);
   assert(A.cols == B.rows);
   int Arows = A.rows;
   int Acols = A.cols;

   for(int j = 0; j < Acols; ++j)
      for(int i = 0; i < Arows; ++i)
         C_ELEM_ID(B, j, i) = C_ELEM_ID(A, i, j);
}

CMatrix CMatrix_Transpose(CMatrix A, trans_type tt)
{
   CMatrix A_TR = CMatrix_init(A.cols, A.rows);
   int Acols = A.cols;
   int Arows = A.rows;

   for(int j = 0; j < Acols; ++j)
      for(int i = 0; i < Arows; ++i)
         A_TR.cid[i][j] = A.cid[j][i];

   if(tt == move) CMatrix_free(A);
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

void CMatrixPrintToFile(CMatrix M, char *info)
{
   FILE *file;
   char str[80];
   sprintf(str, "../results/Matrix_%ix%i_%s.txt", M.rows, M.cols, info);
   file = fopen(str, "w");

   for(int i = 0; i < M.rows; ++i)
   {
      for(int j = 0; j < M.cols; ++j)
         fprintf(file, "%.10e ", C_ELEM_ID(M, i, j));
      fprintf(file, "\n");
   }
   fclose(file);
}
