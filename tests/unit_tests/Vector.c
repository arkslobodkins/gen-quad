/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */

#include "Vector.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <omp.h>

VectorInt VectorInt_init(int n)
{
   assert(n >= 1);

   VectorInt V = {0};
   V.len = n;
   V.id = (int *)calloc(n, sizeof(int));

   return V;
}

void VectorInt_realloc(int n, VectorInt *V)
{
   assert(n >= 1);
   V->id = (int *)realloc(V->id, n*sizeof(int));
   V->len = n;
}

void VectorInt_Assign(VectorInt v1, VectorInt v2)
{
   assert(v1.len == v2.len);
   memcpy(v2.id, v1.id, v1.len*sizeof(int));
}

void VectorInt_free(VectorInt V)
{
   if(V.id != NULL) { free(V.id); V.id = NULL; }
}



Vector Vector_init(int n)
{
   assert(n >= 1);

   Vector V = {0};
   V.len = n;
   V.id = (double *)calloc(n, sizeof(double));
   V.ompId = NULL;

   return V;
}

void Vector_realloc(int n, Vector *V)
{
   assert(n >= 1);
   V->id = (double *)realloc(V->id, n*sizeof(double));
   V->len = n;
}

void Vector_Assign(Vector v1, Vector v2)
{
   assert(v1.len == v2.len);
   memcpy(v2.id, v1.id, v1.len*sizeof(double));
}

void Vector_free(Vector V)
{
   if(V.id != NULL) { free(V.id); V.id = NULL; }
}

#ifdef _OPENMP
void AllocVectorOmpData(Vector *v)
{
   int max_threads = omp_get_max_threads();
   v->ompId = (double **)malloc(max_threads*sizeof(double *));
   for(int i = 0; i < max_threads; ++i)
      v->ompId[i] = (double *)malloc(v->len*sizeof(double));
}

void FreeVectorOmpData(Vector v)
{
   if(v.ompId == NULL) return;
   for(int i = 0; i < omp_get_max_threads(); ++i)
      free(v.ompId[i]);
   free(v.ompId); v.ompId = NULL;
}
#endif



void VPrint(Vector V)
{
   for(int i = 0; i < V.len; ++i)
      printf("V[%i] = %.6e\n", i, V.id[i]);
   printf("\n");
}

double VDot(Vector a, Vector b)
{
   assert(a.len == b.len && b.len > 0);

   double dot = 0.0;
   for(int i = 0; i < a.len; ++i)
      dot += a.id[i]*b.id[i];

   return dot;
}

void VSetToOne(Vector v)
{
   for(int k = 0; k < v.len; ++k)
      v.id[k] = 1.0;
}

void VectorScale(double c, Vector V)
{
   int spacing = 1;
   dscal_(&V.len, &c, V.id, &spacing);
}

void Vector_daxpy(double a, Vector x, Vector y)
{
   assert(x.len == y.len);
   int incx = 1;
   int incy = 1;
   daxpy_(&x.len, &a, x.id, &incx, y.id, &incy);
}

void double_daxpy(int len, double a, double *x, double *y)
{
   int incx = 1;
   int incy = 1;
   daxpy_(&len, &a, x, &incx, y, &incy);
}

void VectorAddScale(double c1, Vector V1, double c2, Vector V2, Vector V3)
{
   assert(V1.len == V2.len && V2.len == V3.len && V3.len > 0);
   for(int i = 0; i < V3.len; ++i)
      V3.id[i] = c1*V1.id[i] + c2*V2.id[i];
}

VMin VectorMin(Vector v)
{
   VMin vMin = {0};

   vMin.min_index = 0;
   vMin.min_value = v.id[0];
   for(int i = 1; i < v.len; ++i)
      if(v.id[i] < vMin.min_value)
      {
         vMin.min_value = v.id[i];
         vMin.min_index = i;
      }

   return vMin;
}

void VectorRemoveElement(int index, Vector *z)
{
   assert(index <= z->len-1);

   for(int i = index; i < z->len-1; ++i)
      z->id[i] = z->id[i+1];

   Vector_realloc(z->len-1, z);
}

double V_ScaledTwoNorm(Vector z)
{
   assert(z.len > 0);
   int incz = 1;
   return dnrm2_(&z.len, z.id, &incz)/sqrt(z.len);
}

double V_TwoNorm(Vector z)
{
   assert(z.len > 0);
   int incz = 1;
   return dnrm2_(&z.len, z.id, &incz);
}

double V_InfNorm(Vector z)
{
   assert(z.len > 0);
   int incz = 1;
   int maxIndex = idamax_(&z.len, z.id, &incz)-1; // subtract 1 (converting from FORTRAN)
   return fabs(z.id[maxIndex]);
}

bool V_CheckNan(Vector z)
{
   assert(z.len > 0);

   for(int i = 0; i < z.len; ++i)
      if(isnan(z.id[i]))
         return true;

   return false;
}

bool V_CheckInf(Vector z)
{
   assert(z.len > 0);

   for(int i = 0; i < z.len; ++i)
      if(isinf(z.id[i]))
         return true;

   return false;
}

void double_dcopy(int n, double *x, double *y)
{
   int incx = 1, incy = 1;
   dcopy_(&n, x, &incx, y, &incy);
}


