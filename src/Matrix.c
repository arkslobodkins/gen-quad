/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */

#include "Matrix.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

RMatrix RMatrix_init(int nRows, int nCols)
{
   RMatrix M = {0};

   M.rows = nRows;
   M.cols = nCols;
   M.len  = nRows*nCols;
   M.id   = (double *)calloc(nRows*nCols, sizeof(double) );
   M.rid  = (double **)malloc(nRows*sizeof(double*) );

   for(int i = 0; i < nRows; ++i)
      M.rid[i] = &M.id[nCols*i];

   return M;
}


void RMatrix_realloc(int nRows, int nCols, RMatrix *M)
{
   M->rid = (double **)realloc(M->rid, nRows*sizeof(double*) );
   M->id  = (double *)realloc(M->id, nRows*nCols*sizeof(double) );
   for(int i = 0; i < nRows; ++i)
      M->rid[i] = &M->id[nCols*i];

   M->rows = nRows;
   M->cols = nCols;
   M->len  = nRows*nCols;
}


void RMatrix_free(RMatrix M)
{
   if(M.id != NULL)  { free(M.id); M.id = NULL; }
   if(M.rid != NULL) { free(M.rid); M.rid = NULL; }
}


void RMatVec(RMatrix M, Vector V, Vector O)
{
   assert(M.cols == V.len);
   assert(O.len >= M.rows);

   int i, j;

   for(i = 0; i < M.rows; ++i)
      for(j = 0; j < M.cols; ++j)
         O.id[i] = 0.0;

   for(i = 0; i < M.rows; ++i)
      for(j = 0; j < M.cols; ++j)
         O.id[i] += R_ELEM_ID(M, i, j) * V.id[j];
}


void PrintRMatrix(const RMatrix M)
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
   M.id   = (double *)calloc(nRows*nCols, sizeof(double) );
   M.cid  = (double **)malloc(nCols*sizeof(double*) );

   for(int j = 0; j < nCols; ++j)
      M.cid[j] = &M.id[nRows*j];
   return M;
}

void CMatrix_Transpose(CMatrix A, CMatrix B)
{
   assert(A.rows == B.cols);
   assert(A.cols == B.rows);

   for(int j = 0; j < A.cols; ++j)
      for(int i = 0; i < A.rows; ++i)
         B.cid[i][j] = A.cid[j][i];
}


void CMatrix_realloc(int nRows, int nCols, CMatrix *M)
{
   M->id  = (double *)realloc(M->id, nRows*nCols*sizeof(double) );
   M->cid = (double **)realloc(M->cid, nCols*sizeof(double*) );
   for(int j = 0; j < nCols; ++j)
      M->cid[j] = &M->id[nRows*j];

   M->rows = nRows;
   M->cols = nCols;
   M->len  = nRows*nCols;
}


void CMatrix_free(CMatrix M)
{
   if(M.id != NULL)  { free(M.id); M.id = NULL; }
   if(M.cid != NULL) { free(M.cid); M.cid = NULL; }
}


void ColMatVec(CMatrix M, Vector V, Vector O)
{
   assert(M.cols == V.len);
   assert(O.len >= M.rows);

   int i, j;

   for(i = 0; i < M.rows; ++i)
      for(j = 0; j < M.cols; ++j)
         O.id[i] = 0.0;

   for(j = 0; j < M.cols; ++j)
      for(i = 0; i < M.rows; ++i)
         O.id[i] += C_ELEM_ID(M, i, j) * V.id[j];

}
