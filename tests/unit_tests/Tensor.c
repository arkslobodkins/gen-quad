#include "Tensor.h"
#include <stdlib.h>

Tensor1D Tensor1D_init(int n)
{
   Tensor1D T = {0};
   T.len = n;
   T.id = (double *)calloc(T.len, sizeof(double));
   return T;
}

Tensor1D VectorToTensor1D(Vector V)
{
   Tensor1D T = {0};
   T.len = V.len;
   T.id = V.id;
   return T;
}

Tensor1D DoubleToTensor1D(int len, double *d)
{
   Tensor1D T = {0};
   T.len = len;
   T.id = d;
   return T;
}

void Tensor1D_free(Tensor1D T)
{
   free(T.id);
}

Tensor2D Tensor2D_init(int dim1, int dim2)
{
   Tensor2D T = {0};
   T.dim1 = dim1;
   T.dim2 = dim2;
   T.len = dim1*dim2;
   T.id = (double *)calloc(T.len, sizeof(double));
   return T;
}

Tensor2D VectorToTensor2D(int dim1, int dim2, Vector V)
{
   Tensor2D T = {0};
   if(dim1*dim2 != V.len)
      return T;

   T.dim1 = dim1;
   T.dim2 = dim2;
   T.len = dim1*dim2;
   T.id = V.id;
   return T;
}

Tensor2D DoubleToTensor2D(int dim1, int dim2, double *d)
{
   Tensor2D T = {0};
   T.dim1 = dim1;
   T.dim2 = dim2;
   T.len = dim1*dim2;
   T.id = d;
   return T;
}

void Tensor2D_free(Tensor2D T)
{
   free(T.id);
}
