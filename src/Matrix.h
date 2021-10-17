#ifndef MATRIX_H
#define MATRIX_H

#include "Vector.h"

#ifdef __cplusplus
extern "C" {
#endif

struct Matrix
{
   int rows;
   int cols;
   double **id;
};
typedef struct Matrix Matrix;

Matrix Matrix_init(int nRows, int nCols);
void Matrix_realloc(int nRows, int nCols, Matrix *M);
void Matrix_free(Matrix M);
void MatVec(Matrix M, Vector V,    Vector O);
void PrintMatrix(const Matrix M);

#ifdef __cplusplus
}
#endif

#endif

