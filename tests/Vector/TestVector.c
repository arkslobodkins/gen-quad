#include "../../src/Vector.h"

#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

const int vlen = 100; // all tests are performed with this length

////////////////////////////////////////////////////////
// Auxiliary routines
////////////////////////////////////////////////////////
bool within_tol(double x, double y)
{
   double tol = pow(10, -15);
   if(fabs(x) < tol && fabs(y) < tol) return true;

   return fabs(x - y) / fmax(x, y) < tol;
}

bool vector_within_tol(Vector a, Vector b)
{
   assert(a.len == b.len);
   for(int i = 0; i < a.len; ++i)
      if( !within_tol(a.id[i], b.id[i]) )
         return false;

   return true;
}

bool vector_equal(Vector a, Vector b)
{
   assert(a.len == b.len);
   for(int i = 0; i < a.len; ++i)
      if(a.id[i] != b.id[i])
         return false;

   return true;
}

Vector initialize_random(double range)
{
   Vector x = Vector_init(vlen);
   for(int i = 0; i < x.len; ++i)
      x.id[i] = range * (double)rand() / (double)RAND_MAX;
   return x;
}

void reinitialize_random(Vector x, double range)
{
   for(int i = 0; i < x.len; ++i)
      x.id[i] = range * (double)rand() / (double)RAND_MAX;
}

////////////////////////////////////////////////////////
// beginning of unit tests
////////////////////////////////////////////////////////
void assert_equal()
{
   Vector a = initialize_random(5.0);

   Vector b = Vector_init(vlen);
   Vector_Assign(a, b);
   assert(vector_equal(a, b));

   Vector_free(a);
   Vector_free(b);

   printf("passed assert_equal\n");
}

void assert_zeros()
{
   Vector x = Vector_init(vlen);

   for(int i = 0; i < x.len; ++i)
      assert(x.id[i] == 0);

   Vector_free(x);

   printf("passed assert_zeros\n");
}

void assert_ones()
{
   Vector x = Vector_init(vlen);
   VSetToOne(x);
   for(int i = 0; i < vlen; ++i)
      assert(x.id[i] == 1);

   Vector_free(x);

   printf("passed assert_ones\n");
}

void assert_dot_itself()
{
   Vector x = initialize_random(5.0);

   double sqsum = 0.0;
   for(int i = 0; i < vlen; ++i)
      sqsum += x.id[i] * x.id[i];

   double vdot = VDot(x, x);
   assert(within_tol(vdot, sqsum));

   Vector_free(x);

   printf("passed assert_dot_itself\n");
}

void assert_scale()
{
   Vector a = initialize_random(5.0);
   Vector b = Vector_init(vlen);
   Vector_Assign(a, b);

   VScale(1.0, b);
   vector_within_tol(a, b);

   Vector_Assign(a, b);
   VScale(-1.0, b);
   for(int i = 0; i < a.len; ++i)
      assert( within_tol(-a.id[i], b.id[i]) );

   Vector_Assign(a, b);
   VScale(2.0, b);
   for(int i = 0; i < a.len; ++i)
      assert( within_tol(a.id[i] + a.id[i], b.id[i]) );

   Vector_free(a);
   Vector_free(b);

   printf("passed assert_scale\n");
}

void assert_daxpy()
{
   Vector a = initialize_random(5.0);
   Vector b = Vector_init(vlen);
   Vector_Assign(a, b);

   Vector_Assign(a, b);           // b = a
   Vector_daxpy(2.0, a, b);       // b = 2*a + a = 3*a
   for(int i = 0; i < vlen; ++i)
      assert( within_tol(a.id[i] + a.id[i] + a.id[i], b.id[i]) );

   Vector_free(a);
   Vector_free(b);

   printf("passed assert_daxpy\n");
}

void assert_add_scale()
{
   Vector a = initialize_random(5.0);
   Vector b = initialize_random(5.0);
   Vector c = Vector_init(vlen);

   VAddScale(1.0, a, 1.0, b, c);
   for(int i = 0; i < vlen; ++i)
      assert( within_tol(a.id[i] + b.id[i], c.id[i]) );

   reinitialize_random(a, 5.0);
   reinitialize_random(b, 5.0);

   VAddScale(1.0, a, 2.0, b, c);
   for(int i = 0; i < vlen; ++i)
      assert( within_tol(a.id[i] + b.id[i] + b.id[i], c.id[i]) );

   reinitialize_random(a, 5.0);
   VAddScale(1.0, a, -1.0, a, c);
   for(int i = 0; i < vlen; ++i)
      assert(within_tol(0.0, c.id[i]));

   Vector_free(a);
   Vector_free(b);
   Vector_free(c);

   printf("passed assert_add_scale\n");
}

void assert_min()
{
   Vector x = Vector_init(vlen);

   for(int i = 0; i < x.len; ++i)
      x.id[i] = i;
   x.id[x.len-1] = -10.;

   VMin min = VectorMin(x);
   assert(min.min_index == x.len-1);
   assert(min.min_value == x.id[x.len-1]);

   x.id[0] = x.id[x.len-1]-0.0001;
   min = VectorMin(x);
   assert(min.min_index == 0);
   assert(min.min_value == x.id[0]);

   Vector_free(x);

   printf("passed assert_min\n");
}

void assert_max()
{
   Vector x = Vector_init(vlen);

   for(int i = 0; i < x.len; ++i)
      x.id[i] = i;
   x.id[x.len-1] = x.len*3.0;

   VMax max = VectorMax(x);
   assert(max.max_index == x.len-1);
   assert(max.max_value == x.id[x.len-1]);

   x.id[0] = x.id[x.len-1]+0.0001;
   max = VectorMax(x);
   assert(max.max_index == 0);
   assert(max.max_value == x.id[0]);

   Vector_free(x);

   printf("passed assert_max\n");
}

void assert_remove()
{
   Vector a = Vector_init(vlen);
   Vector b = Vector_init(vlen-1);

   for(int i = 0; i < a.len; ++i) a.id[i] = i;
   for(int i = 0; i < b.len; ++i) b.id[i] = i;
   VRemoveElement(a.len-1, &a);
   assert(vector_equal(a, b));

   Vector_realloc(vlen, &a);
   for(int i = 0; i < a.len; ++i) a.id[i] = i-1;
   VRemoveElement(0, &a);
   assert(vector_equal(a, b));

   Vector_realloc(vlen, &a);
   for(int i = 0; i <= a.len/2; ++i) a.id[i] = i;
   for(int i = a.len/2+1; i < a.len; ++i) a.id[i] = i-1;
   VRemoveElement(a.len/2, &a);
   assert(vector_equal(a, b));

   Vector_free(a);
   Vector_free(b);

   printf("passed assert_remove\n");
}

void assert_increasing()
{
   Vector x = Vector_init(vlen);
   for(int i = 0; i < x.len; ++i) x.id[i] = i;
   assert(V_IsIncreasing(x));

   x.id[0] = 5;
   assert(!V_IsIncreasing(x));

   Vector_free(x);
   printf("passed assert_increasing\n");
}

int main()
{
   assert_zeros();
   assert_equal();
   assert_dot_itself();
   assert_ones();
   assert_scale();
   assert_daxpy();
   assert_add_scale();
   assert_min();
   assert_max();
   assert_remove();
   assert_increasing();

   return EXIT_SUCCESS;
}
