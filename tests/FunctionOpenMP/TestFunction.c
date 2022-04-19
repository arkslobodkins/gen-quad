#include "../../src/GetFunction.h"
#include "../../src/Quadrature.h"
#include "../../src/Vector.h"
#include "../../src/get_time.h"

#include <stdlib.h>
#include <stdio.h>

#include <mkl.h>
#include <omp.h>

int main()
{
   mkl_set_threading_layer(MKL_THREADING_SEQUENTIAL);
   printf("testing with %i threads\n", omp_get_max_threads());

   int numTests = 10;
   // set problem parameters, typical 4-d size problem
   int nodes = 1200;
   int dim = 4;
   int dims[1] = {1};
   dims[0] = dim;
   int deg = 15;
   DOMAIN_TYPE D = CUBE;

   // initialize elements to random values in [0, 1], similar to actual quadrature values
   quadrature *q[numTests];
   for(int i = 0; i < numTests; ++i)
   {
      q[i] = quadrature_init_full(nodes, dim, dims, deg, D);
      quadrature_fill_random(q[i]);
   }

   Vector fSerial = Vector_init(q[0]->basis->numFuncs);
   Vector fParallel = Vector_init(q[0]->basis->numFuncs);
   AllocVectorOmpData(&fParallel, omp_get_max_threads());
   for(int i = 0; i < numTests; ++i)
      QuadAllocBasisOmp(q[i], omp_get_max_threads());

   ////////////////////////////////////////////////////////////
   // serial and parallel norm and times
   printf("\n");
   for(int i = 0; i < numTests; ++i)
   {
      TIME(GetFunction(q[i], fSerial));
      TIME(GetFunctionOmp(q[i], fParallel));
      printf("maximum absolute difference = %.16e\n", VectorMaxDifference(fSerial, fParallel));
      printf("maximum relative difference = %.16e\n", VectorMaxRelativeDifference(fSerial, fParallel));
      printf("\n");
   }
   printf("\n");
   ////////////////////////////////////////////////////////////

   Vector_free(fSerial);
   FreeVectorOmpData(fParallel);
   Vector_free(fParallel);

   for(int i = 0; i < numTests; ++i)
      QuadFreeBasisOmp(q[i]);
   for(int i = 0; i < numTests; ++i)
      quadrature_free(q[i]);

   return EXIT_SUCCESS;
}
