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
void assert_equal(int line)
{
   Vector a = initialize_random(5.0);

   Vector b = Vector_init(vlen);
   Vector_Assign(a, b);
   assert(vector_equal(a, b));

   Vector_free(a);
   Vector_free(b);

   printf("passed assert_equal on line %i\n", line);
}

void assert_zeros(int line)
{
   Vector x = Vector_init(vlen);

   for(int i = 0; i < x.len; ++i)
      assert(x.id[i] == 0);

   Vector_free(x);

   printf("passed assert_zeros on line %i\n", line);
}

void assert_ones(int line)
{
   Vector x = Vector_init(vlen);
   VSetToOne(x);
   for(int i = 0; i < vlen; ++i)
      assert(x.id[i] == 1);

   Vector_free(x);

   printf("passed assert_ones on line %i\n", line);
}

void assert_dot_itself(int line)
{
   Vector x = initialize_random(5.0);

   double sqsum = 0.0;
   for(int i = 0; i < vlen; ++i)
      sqsum += x.id[i] * x.id[i];

   double vdot = VDot(x, x);
   assert(within_tol(vdot, sqsum));

   Vector_free(x);

   printf("passed assert_dot_itself on line %i\n", line);
}

void assert_scale(int line)
{
   Vector a = initialize_random(5.0);
   Vector b = Vector_init(vlen);
   Vector_Assign(a, b);

   VScale(1.0, b);
   vector_within_tol(a, b);

   Vector_Assign(a, b);
   VScale(-1.0, b);
   for(int i = 0; i < a.len; ++i)
      assert(within_tol(-a.id[i], b.id[i]));

   Vector_Assign(a, b);
   VScale(2.0, b);
   for(int i = 0; i < a.len; ++i)
      assert(within_tol(a.id[i] + a.id[i], b.id[i]));

   Vector_free(a);
   Vector_free(b);

   printf("passed assert_scale on line %i\n", line);
}

void assert_daxpy(int line)
{
   Vector a = initialize_random(5.0);
   Vector b = Vector_init(vlen);
   Vector_Assign(a, b);

   Vector_Assign(a, b);           // b = a
   Vector_daxpy(2.0, a, b);       // b = 2*a + a = 3*a
   for(int i = 0; i < vlen; ++i)
      assert(within_tol(a.id[i] + a.id[i] + a.id[i], b.id[i]));

   Vector_free(a);
   Vector_free(b);

   printf("passed assert_daxpy on line %i\n", line);
}

void assert_add_scale(int line)
{
   Vector a = initialize_random(5.0);
   Vector b = initialize_random(5.0);
   Vector c = Vector_init(vlen);

   VAddScale(1.0, a, 1.0, b, c);
   for(int i = 0; i < vlen; ++i)
      assert(within_tol(a.id[i] + b.id[i], c.id[i]));

   reinitialize_random(a, 5.0);
   reinitialize_random(b, 5.0);

   VAddScale(1.0, a, 2.0, b, c);
   for(int i = 0; i < vlen; ++i)
      assert(within_tol(a.id[i] + b.id[i] + b.id[i], c.id[i]));

   reinitialize_random(a, 5.0);
   VAddScale(1.0, a, -1.0, a, c);
   for(int i = 0; i < vlen; ++i)
      assert(within_tol(0.0, c.id[i]));

   Vector_free(a);
   Vector_free(b);
   Vector_free(c);

   printf("passed assert_add_scale on line %i\n", line);
}

int main(int argc, char *argv[])
{
   assert_zeros(__LINE__);
   assert_equal(__LINE__);
   assert_dot_itself(__LINE__);
   assert_ones(__LINE__);
   assert_scale(__LINE__);
   assert_daxpy(__LINE__);
   assert_add_scale(__LINE__);

   return EXIT_SUCCESS;
}
