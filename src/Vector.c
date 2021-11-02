/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */

#include "Vector.h"
#include "GENERAL_QUADRATURE.h"

#include <assert.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <math.h>


Vector_bool Vector_bool_init(int n)
{
   Vector_bool V = {0};

   V.len = n;
   V.id = (bool *)calloc(n, sizeof(bool));

   return V;
}

void Vector_bool_realloc(int n, Vector_bool *V)
{
   V->len = n;
   V->id = (bool *)realloc(V->id, n*sizeof(bool));
}

void Vector_bool_free(Vector_bool V)
{
   if(V.id != NULL) {
      free(V.id);
      V.id = NULL;
   }
}


Vector_int Vector_int_init(int n)
{
   Vector_int V = {0};
   V.len = n;
   V.id = (int *)calloc(n, sizeof(int));

   return V;
}

void Vector_int_realloc(int n, Vector_int *V)
{
   V->len = n;
   V->id = (int *)realloc(V->id, n*sizeof(int));
}

void Vector_int_free(Vector_int V)
{
   if(V.id != NULL) {
      free(V.id);
      V.id = NULL;
   }
}


Vector Vector_init(int n)
{
   Vector V = {0};
   V.len = n;
   V.id = (double *)calloc(n, sizeof(double));

   return V;
}

void Vector_realloc(int n, Vector *V)
{
   V->len = n;
   V->id = (double *)realloc(V->id, n*sizeof(double));
}

void Vector_Assign(const Vector v1, Vector v2)
{
   assert(v1.len == v2.len);
   for(int i = 0; i < v1.len; ++i)
      v2.id[i] = v1.id[i];
}

void Vector_free(Vector V)
{
   if(V.id != NULL) {
      free(V.id);
      V.id = NULL;
   }
}

void Vector_Print(const Vector V)
{
   for(int i = 0; i < V.len; ++i)
      printf("V[%i] = %.6e\n", i, V.id[i]);
   printf("\n");
}

double VDot(const Vector a, const Vector b)
{
   assert(a.len == b.len && b.len > 0);

   double dot = 0.0;
   for(int i = 0; i < a.len; ++i)
      dot += a.id[i]*b.id[i];

   return dot;
}

void VectorAddScale(double c1, const Vector V1, double c2, const Vector V2, Vector V3)
{
   assert(V1.len == V2.len  && V2.len == V3.len  && V3.len > 0);
   for(int i = 0; i < V3.len; ++i)
      V3.id[i] = c1*V1.id[i] + c2*V2.id[i];
}

VMin VectorMin(const Vector v)
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


void VRemoveElement(int index, Vector *z)
{
   assert(index <= z->len-1);

   for(int i = index; i < z->len-1; ++i)
      z->id[index] = z->id[i+1];

   Vector_realloc(z->len-1, z);
}

void VIntRemoveElement(int index, Vector_int *z)
{
   assert(index <= z->len-1);

   for(int i = index; i < z->len-1; ++i)
      z->id[index] = z->id[i+1];

   Vector_int_realloc(z->len-1, z);
}


double V_ScaledTwoNorm(const Vector z)
{
   assert(z.len > 0);

   double norm = 0.0;
   for(int i = 0; i < z.len; ++i)
      norm += SQUARE(z.id[i]);

   norm = SQRT(norm/z.len);
   return norm;
}

double V_TwoNorm(const Vector z)
{
   assert(z.len > 0);

   double norm = 0.0;
   for(int i = 0; i < z.len; ++i)
      norm += SQUARE(z.id[i]);

   norm = SQRT(norm);
   return norm;
}

double V_InfNorm(const Vector z)
{
   assert(z.len > 0);

   double norm = 0.0;
   for(int i = 0; i < z.len; ++i)
      norm = MAX(fabs(z.id[i]), norm);

   return norm;
}

bool V_CheckInfAndNan(const Vector z)
{
   assert(z.len > 0);

   for(int i = 0; i < z.len; ++i)
      if( (isnan(z.id[i]) == 1) || (isinf(z.id[i]) == 1) )
         return true;

   return false;
}

int IntPower(int x, int power)
{
   int result = 1;
   for(int i = 0; i < power; ++i)
      result *= x;

   return result;
}
