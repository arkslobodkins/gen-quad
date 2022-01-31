#ifndef MATRIX_H
#define MATRIX_H

#include "Vector.h"

#ifdef __cplusplus
extern "C" {
#endif

#define R_ELEM_ID(A, i, j) ( (A.rid)[i][j] )
#define C_ELEM_ID(A, i, j) ( (A.cid)[j][i] )

typedef struct
{
   int rows;
   int cols;
   int len;
   double *id;
   double **rid;
} RMatrix;

typedef struct CMatrix
{
   int rows;
   int cols;
   int len;
   double *id;
   double **cid;
} CMatrix;

RMatrix RMatrix_init(int nRows, int nCols);
void RMatrix_realloc(int nRows, int nCols, RMatrix *M);
void RMatrix_free(RMatrix M);
void RMatVec(RMatrix M, Vector x, Vector y);
void RMatrixPrint(RMatrix M);

CMatrix CMatrix_init(int nRows, int nCols);
void CMatrix_realloc(int nRows, int nCols, CMatrix *M);
void CMatrix_free(CMatrix M);
void CMatrix_Assign(CMatrix A, CMatrix B);
void CMatrix_LoadColumn(int col, Vector x, CMatrix M);
void CMatrix_LoadColumnDD(int col, int len, double *x, CMatrix M);
void CMatrix_GetRow(int row, CMatrix M, Vector x);
void CMatrix_GetColumn(int col, CMatrix M, Vector x);
void CMatrix_ColumnScale(int col, double c, CMatrix M);
CMatrix CMatrix_Transpose(CMatrix A);
void CMatVec(CMatrix M, Vector x, Vector y);
void CMatrixPrint(CMatrix M);
void CMatrixPrintToFile(CMatrix M, char *info);

#ifdef __cplusplus
}
#endif

#endif

