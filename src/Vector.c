/* Arkadijs Slobodkins
 * SMU Mathematics
 * Copyright 2022
 */


#include "Vector.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>

#include <omp.h>
#include <mkl_blas.h>

static bool comparedouble(double x, double y);

VectorInt VectorInt_init(int n)
{
   assert(n >= 1);

   VectorInt z = {0, 0};
   z.len = n;
   z.id = (int *)calloc(n, sizeof(int));

   return z;
}

void VectorInt_realloc(int n, VectorInt *z)
{
   assert(n >= 1);
   z->id = (int *)realloc(z->id, n*sizeof(int));
   z->len = n;
}

void VectorInt_Assign(VectorInt V1, VectorInt V2)
{
   assert(V1.len == V2.len);
   memcpy(V2.id, V1.id, V1.len*sizeof(int));
}

void VectorInt_free(VectorInt z)
{
   if(z.id != NULL) { free(z.id); }
}



Vector Vector_init(int n)
{
   assert(n >= 1);

   Vector z;
   z.len = n;
   z.id = (double *)calloc(n, sizeof(double));
   z.ompId = NULL;
   z.ompLen = 0;

   return z;
}

Vector Vector_uninitialized(int n)
{
   assert(n >= 1);

   Vector z;
   z.len = n;
   z.id = (double *)malloc(n*sizeof(double));
   for(int i = 0; i < n; ++i)
      z.id[i] = -100.;
   z.ompId = NULL;
   z.ompLen = 0;

   return z;
}

void Vector_realloc(int n, Vector *z)
{
   assert(n >= 1);
   z->id = (double *)realloc(z->id, n*sizeof(double));
   z->len = n;
}

void Vector_Assign(Vector V1, Vector V2)
{
   assert(V1.len >= 1);
   assert(V1.len == V2.len);
   memcpy(V2.id, V1.id, V1.len*sizeof(double));
}

void Vector_free(Vector z)
{
   if(z.id != NULL) { free(z.id); }
}

#ifdef _OPENMP
void AllocVectorOmpData(Vector *z, int num_threads)
{
   z->ompLen = num_threads;
   z->ompId = (double **)malloc(num_threads*sizeof(double *));
   for(int i = 0; i < num_threads; ++i)
      z->ompId[i] = (double *)malloc(z->len*sizeof(double));
}

void FreeVectorOmpData(Vector z)
{
   if(z.ompId == NULL) return;
   for(int i = 0; i < z.ompLen; ++i)
      free(z.ompId[i]);
   free(z.ompId);
}
#endif

void VPrint(Vector V)
{
   for(int i = 0; i < V.len; ++i)
      printf("V[%i] = %.6e\n", i, V.id[i]);
   printf("\n");
}

double VDot(Vector x, Vector y)
{
   assert(x.len == y.len && y.len > 0);

   double dot = 0.0;
   for(int i = 0; i < x.len; ++i)
      dot += x.id[i] * y.id[i];

   return dot;
}

void VSetToZero(Vector z)
{
   assert(z.len >= 1);
   memset(z.id, 0, z.len*sizeof(double));
}

void VSetToOne(Vector z)
{
   assert(z.len >= 1);
   for(int k = 0; k < z.len; ++k)
      z.id[k] = 1.0;
}

void VScale(double c, Vector z)
{
   assert(z.len >= 1);
   int spacing = 1;
   dscal(&z.len, &c, z.id, &spacing);
}

void Vector_daxpy(double a, Vector x, Vector y)
{
   assert(x.len == y.len);
   int incx = 1;
   int incy = 1;
   daxpy(&x.len, &a, x.id, &incx, y.id, &incy);
}

void VAdd(Vector V1, Vector V2, Vector V3)
{
   assert(V1.len == V2.len && V2.len == V3.len && V3.len > 0);
   for(int i = 0; i < V3.len; ++i)
      V3.id[i] = V1.id[i] + V2.id[i];
}

void VAddScale(double c1, Vector V1, double c2, Vector V2, Vector V3)
{
   assert(V1.len == V2.len && V2.len == V3.len && V3.len > 0);
   for(int i = 0; i < V3.len; ++i)
      V3.id[i] = c1*V1.id[i] + c2*V2.id[i];
}

VMin VectorMin(Vector z)
{
   assert(z.len > 0);
   VMin vMin = {0, 0};

   vMin.min_index = 0;
   vMin.min_value = z.id[0];
   for(int i = 1; i < z.len; ++i)
      if(z.id[i] < vMin.min_value)
      {
         vMin.min_value = z.id[i];
         vMin.min_index = i;
      }

   return vMin;
}

void VRemoveElement(int index, Vector *z)
{
   assert(index >= 0 && index <= z->len-1);

   for(int i = index; i < z->len-1; ++i)
      z->id[i] = z->id[i+1];

   Vector_realloc(z->len-1, z);
}

VMax VectorMax(Vector z)
{
   assert(z.len > 0);
   VMax vMax;

   vMax.max_index = 0;
   vMax.max_value = z.id[0];
   for(int i = 1; i < z.len; ++i)
      if(z.id[i] > vMax.max_value)
      {
         vMax.max_value = z.id[i];
         vMax.max_index = i;
      }

   return vMax;
}

double VectorMaxDifference(Vector x, Vector y)
{
   assert(x.len == y.len);
   double maxDiff = fabs(x.id[0] - y.id[0]);

   for(int i = 1; i < x.len; ++i)
   {
      double iDiff = fabs(x.id[i] - y.id[i]);
      maxDiff = iDiff > maxDiff ? iDiff : maxDiff;
   }
   return maxDiff;
}


double V_ScaledTwoNorm(Vector z)
{
   assert(z.len > 0);
   int incz = 1;
   return dnrm2(&z.len, z.id, &incz)/sqrt(z.len);
}

double V_TwoNorm(Vector z)
{
   assert(z.len > 0);
   int incz = 1;
   return dnrm2(&z.len, z.id, &incz);
}

double V_InfNorm(Vector z)
{
   assert(z.len > 0);
   int incz = 1;
   int maxIndex = idamax(&z.len, z.id, &incz)-1; // subtract 1 (converting from FORTRAN)
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

bool V_CheckInfOrNan(Vector z)
{
   return V_CheckInf(z) || V_CheckNan(z);
}

bool V_IsUninitialized(Vector z)
{
   for(int i = 0; i < z.len; ++i)
      if(comparedouble(z.id[i], -100.))
         return true;
   return false;
}

bool V_IsIncreasing(Vector v)
{
   for(int i = 1; i < v.len; ++i)
      if(v.id[i] < v.id[i-1])
         return false;

   return true;
}

void double_daxpy(int len, double a, double *x, double *y)
{
   int incx = 1;
   int incy = 1;
   daxpy(&len, &a, x, &incx, y, &incy);
}

void double_dcopy(int n, double *x, double *y)
{
   int incx = 1, incy = 1;
   dcopy(&n, x, &incx, y, &incy);
}

static bool comparedouble(double x, double y)
{
   if( fabs(x-y) < pow(10, -15) ) return true;
   return false;
}


