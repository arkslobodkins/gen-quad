#ifndef TENSOR_H
#define TENSOR_H

#include "Vector.h"

#ifdef __cplusplus
extern "C" {
#endif

#define TID1(T, i) ( T.id[(i)] )
#define TID2(T, i, j) ( T.id[(i)*(T.dim2)+(j)] )

typedef struct
{
   int len;
   double *id;
} Tensor1D;

typedef struct
{
   int dim1;
   int dim2;
   int len;
   double *id;
} Tensor2D;

Tensor1D Tensor1D_init(int n);
Tensor1D VectorToTensor1D(Vector V);
Tensor1D DoubleToTensor1D(int len, double *d);
void Tensor1D_free(Tensor1D T);

Tensor2D Tensor2D_init(int dim1, int dim2);
Tensor2D VectorToTensor2D(int dim1, int dim2, Vector V);
Tensor2D DoubleToTensor2D(int dim1, int dim2, double *d);
void Tensor2D_free(Tensor2D T);

#ifdef __cplusplus
}
#endif

#endif
