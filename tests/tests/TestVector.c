#include "../../src/Vector.h"

#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
   Vector b = Vector_init(8);
   Vector v = Vector_init(8);
   Vector z = Vector_init(8);
   v.id[0] = 1.; v.id[1] = -1.;
   v.id[2] = 2.; v.id[3] = -2.;
   v.id[4] = 3.; v.id[5] = -3.;
   v.id[6] = 4.; v.id[7] = -4.;
   Vector_Assign(v, z);
   Vector_Assign(v, b);
   for(int i = 0; i < v.len; ++i)
      assert(v.id[i] == z.id[i]);


   printf("testing Vector class, initial values:\n");
   VPrint(v);
   printf("VDot(v, v) = %lf\n\n", VDot(v, v));


   VSetToOne(v);
   printf("after setting to one:\n");
   VPrint(v);


   v.id[0] = 1.; v.id[1] = -1.;
   v.id[2] = 2.; v.id[3] = -2.;
   v.id[4] = 3.; v.id[5] = -3.;
   v.id[6] = 4.; v.id[7] = -4.;
   printf("reinitializing to original values:\n");
   VPrint(v);


   VScale(1.5, v);
   printf("scaling by 1.5:\n");
   VPrint(v);


   Vector_daxpy(2.0, v, z);
   for(int i = 0; i < v.len; ++i)
      assert(fabs(z.id[i] - 4.0*b.id[i]) < pow(10, -15));


   VAddScale(-2.0, v, 2.0, b, z);
   for(int i = 0; i < v.len; ++i)
      assert( z.id[i] - (-2.*v.id[i] + 2.*b.id[i]) < pow(10, -15) );


   VRemoveElement(3, &v);
   printf("removing element #3:\n");
   printf("number of elements = %i\n", v.len);
   VPrint(v);


   VMax vmax = VectorMax(v);
   printf("max index = %i\n", vmax.max_index);
   printf("max value = %lf\n", vmax.max_value);


   VMin vmin = VectorMin(v);
   printf("min index = %i\n", vmin.min_index);
   printf("min value = %lf\n\n", vmin.min_value);


   printf("scaled two norm of v = %lf\n", V_ScaledTwoNorm(v));
   printf("two norm of v = %lf\n", V_TwoNorm(v));
   printf("infinity norm of v = %lf\n\n", V_InfNorm(v));


   v.id[1] = pow(10., 1000);
   VPrint(v);
   printf("vector isinf after setting an entry to infinity: %d \n", V_CheckInf(v));
   printf("vector isinf or isnan after setting an entry to infinity: %d \n", V_CheckInfOrNan(v));
   VPrint(v);
   v.id[1] = NAN;
   printf("vector isnan after setting an entry to nan: %d \n", V_CheckNan(v));
   printf("vector isinf or isnan after setting an entry to nan: %d \n", V_CheckInfOrNan(v));
   VPrint(v);
   v.id[1] = -100.;
   printf("checking for uninitialized value: %d \n", V_IsUninitialized(v));
   VPrint(v);


   double a1[3] = {1.0, 2.0, 3.0};
   double a2[3] = {2.0, 3.0, 4.0};
   double_daxpy(3, 2.0, a1, a2);
   for(int i = 0; i < 3; ++i)
      assert(a2[i] == 3.0*a1[i]+1.0);


   double_dcopy(v.len, v.id, z.id);
   for(int i = 0; i < v.len; ++i)
      assert(v.id[i] == z.id[i]);


   Vector_free(b);
   Vector_free(v);
   Vector_free(z);
   return EXIT_SUCCESS;
}
