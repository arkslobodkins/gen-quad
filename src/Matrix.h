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
void PrintRMatrix(const RMatrix M);

CMatrix CMatrix_init(int nRows, int nCols);
void CMatrix_realloc(int nRows, int nCols, CMatrix *M);
void CMatrix_free(CMatrix M);
void CMatrix_Assign(const CMatrix A, CMatrix B);
void CMatrix_Transpose(const CMatrix A, CMatrix B);
void CMatVec(const CMatrix M, const Vector x, Vector y);

#ifdef __cplusplus
}
#endif

#endif

