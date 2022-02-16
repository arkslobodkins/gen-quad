#include "Tensor.h"
#include <stdlib.h>

Tensor Tensor_init(int dim1, int dim2)
{
   Tensor T = {0};
   T.dim1 = dim1;
   T.dim2 = dim2;
   T.len = dim1*dim2;
   T.id = (double *)calloc(T.len, sizeof(double));
   return T;
}

Tensor VectorToTensor(int dim1, int dim2, Vector V)
{
   Tensor T = {0};
   if(dim1*dim2 != V.len)
      return T;

   T.dim1 = dim1;
   T.dim2 = dim2;
   T.len = dim1*dim2;
   T.id = V.id;
   return T;
}

Tensor DoubleToTensor(int dim1, int dim2, double *d)
{
   Tensor T = {0};
   T.dim1 = dim1;
   T.dim2 = dim2;
   T.len = dim1*dim2;
   T.id = d;
   return T;
}

void Tensor_free(Tensor T)
{
   free(T.id);
}
