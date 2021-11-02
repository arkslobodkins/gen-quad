#include <stdlib.h>
#include <gtest/gtest.h>
#include "../src/GENERAL_QUADRATURE.h"

int main(int argc, char *argv[])
{
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}


TEST(GEN_QUAD_TEST, GENERAL_QUADRATURE_MACROS)
{
   ASSERT_EQ(ONE, 1);
   ASSERT_EQ(ONE+1, 2);
   ASSERT_EQ(6/(ONE+1), 3);

   ASSERT_EQ(POW_INT(1, 5), 1);
   ASSERT_EQ(POW_INT(1+2, 3-1), 9);
   ASSERT_EQ(POW_INT(1+2, 0), 1);

   ASSERT_EQ(POW_DOUBLE(2.0, 3.0), 8.0);

   ASSERT_EQ(MAX(-1,-1), -1);
   ASSERT_EQ(MAX(-1,2), 2);
   ASSERT_EQ(MAX(-1.2,-200.0), -1.2);
   ASSERT_EQ(MAX(101,303), 303);
   ASSERT_EQ(MAX(5+5, 9), 10);

   ASSERT_EQ(MIN(-1,-1), -1);
   ASSERT_EQ(MIN(-1,2), -1);
   ASSERT_EQ(MIN(-1.2,-200.0), -200.0);
   ASSERT_EQ(MIN(101,303), 101);
   ASSERT_EQ(MIN(5+5, 9), 9);

   ASSERT_EQ(SQUARE(3+5), 64);
   ASSERT_EQ(SQUARE(-2-1), 9);
   ASSERT_EQ(SQUARE((-2-1)/0.5), 36.0);

   ASSERT_EQ(SQRT(4.0), 2.0);

   ASSERT_EQ(QUAD_TOL, pow(10, -15));

   ASSERT_EQ(size_int, sizeof(int));

   ASSERT_EQ(size_double, sizeof(double));

   ASSERT_EQ(SIZE_INT(5), 5*sizeof(int));
   ASSERT_EQ(SIZE_INT(3*(4+2-1)), 15*sizeof(int));

   ASSERT_EQ(SIZE_DOUBLE(5.0), 5.0*sizeof(double));
   ASSERT_EQ(SIZE_DOUBLE(3.0*(4.0+2.0-1.0)), 15.0*sizeof(double));
}


TEST(GEN_QUAD_TEST, VECTOR)
{
   Vector v = Vector_init(3);
   ASSERT_EQ(v.len, 3);
   for(int i = 0; i < v.len; ++i)
      v.id[i] = i;
   for(int i = 0; i < v.len; ++i)
      EXPECT_EQ(v.id[i], i);

   Vector_realloc(5, &v);
   ASSERT_EQ(v.len, 5);
   for(int i = 0; i < 3; ++i)
      EXPECT_EQ(v.id[i], i);

   Vector w = Vector_init(v.len);
   Vector_Assign(v, w);
   for(int i = 0; i < v.len; ++i)
      EXPECT_EQ(v.id[i], w.id[i]);
}


