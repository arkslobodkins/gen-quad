/* Arkadijs Slobodkins
 * SMU Mathematics
 * Copyright 2022
 */


#ifndef TENSOR_H
#define TENSOR_H

#include "Vector.h"
#include <assert.h>

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


static inline void set_tid1(Tensor1D T1, int i, double d)
{
   assert(i >= 0 && i <= T1.len-1);
   T1.id[i] = d;
}

static inline double get_tid1(Tensor1D T1, int i)
{
   assert(i >= 0 && i <= T1.len-1);
   return T1.id[i];
}

static inline void set_tid2(Tensor2D T2, int i, int j, double d)
{
   assert(i >= 0 && i <= T2.dim1-1);
   assert(j >= 0 && j <= T2.dim2-1);
   T2.id[i*T2.dim2+j] = d;
}

static inline double get_tid2(Tensor2D T2, int i, int j)
{
   assert(i >= 0 && i <= T2.dim1-1);
   assert(j >= 0 && j <= T2.dim2-1);
   return T2.id[i*T2.dim2+j];
}


#ifdef __cplusplus
}
#endif

#endif
