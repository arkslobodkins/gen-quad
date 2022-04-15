#include "../../src/GetFunction.h"
#include "../../src/Quadrature.h"
#include "../../src/Vector.h"
#include "../../src/get_time.h"

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

int main()
{
   int numTests = 5;
   // set problem parameters, typical 4-d size problem
   int nodes = 700;
   int dim = 4;
   int dims[1] = {1};
   dims[0] = dim;
   int deg = 13;
   DOMAIN_TYPE D = CUBE;

   // initialize elements to random values in [0, 1], similar to actual quadrature values
   quadrature *q[numTests];
   for(int i = 0; i < numTests; ++i)
   {
      q[i] = quadrature_init_full(nodes, dim, dims, deg, D);
      quadrature_fill_random(q[i]);
   }


   ////////////////////////////////////////////////////////////
   // serial norm and times
   Vector fSerial = Vector_init(q[0]->basis->numFuncs);
   printf("\n");
   printf("timing serial GetFunction\n");
   for(int i = 0; i < numTests; ++i)
   {
      TIME(GetFunction(q[i], fSerial));
      printf("TwoNorm of GetFunction = %.16e\n", V_TwoNorm(fSerial));
      printf("\n");
   }
   printf("\n");
   Vector_free(fSerial);
   ////////////////////////////////////////////////////////////


   ////////////////////////////////////////////////////////////
   // parallel norm and times
   Vector fParallel = Vector_init(q[0]->basis->numFuncs);
   AllocVectorOmpData(&fParallel, omp_get_max_threads());
   for(int i = 0; i < numTests; ++i)
      QuadAllocBasisOmp(q[i], omp_get_max_threads());

   printf("timing parallel GetFunctionOmp\n");
   for(int i = 0; i < numTests; ++i)
   {
      TIME(GetFunctionOmp(q[i], fParallel));
      printf("TwoNorm of GetFunctionOmp = %.16e\n", V_TwoNorm(fParallel));
      printf("\n");
   }
   printf("\n");
   FreeVectorOmpData(fParallel);
   Vector_free(fParallel);
   ////////////////////////////////////////////////////////////


   for(int i = 0; i < numTests; ++i)
      QuadFreeBasisOmp(q[i]);
   for(int i = 0; i < numTests; ++i)
      quadrature_free(q[i]);

   return EXIT_SUCCESS;
}
