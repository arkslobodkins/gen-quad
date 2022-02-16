#ifndef TENSOR_H
#define TENSOR_H

#include "Vector.h"

#ifdef __cplusplus
extern "C" {
#endif

#define TID(T, i, j) ( T.id[(i)*(T.dim2)+(j)] )

typedef struct
{
   int dim1;
   int dim2;
   int len;
   double **tid;
   double *id;
   void(*elem)(int, int);
} Tensor;

Tensor Tensor_init(int dim1, int dim2);
Tensor VectorToTensor(int dim1, int dim2, Vector V);
Tensor DoubleToTensor(int dim1, int dim2, double *d);
void Tensor_free(Tensor T);

#ifdef __cplusplus
}
#endif

#endif
