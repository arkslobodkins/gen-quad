#include "../../src/GetJacobian.h"
#include "../../src/Quadrature.h"
#include "../../src/Matrix.h"
#include "../../src/get_time.h"
#include "../../src/get_time.h"

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

int main()
{
   int numTests = 5;
   // set problem parameters, typical 4-d size problem
   int nodes = 500;
   int dim = 4;
   int dims[1] = {1};
   dims[0] = dim;
   int deg = 12;
   DOMAIN_TYPE D = CUBE;

   // initialize elements to random values in [0, 1], similar to actual quadrature values
   quadrature *q[numTests];
   for(int i = 0; i < numTests; ++i)
   {
      q[i] = quadrature_init_full(nodes, dim, dims, deg, D);
      quadrature_fill_random(q[i]);
      for(int j = 0; j < q[i]->z.len; ++j)
      {
         q[i]->z.id[j] *= 0.5;
         q[i]->z.id[j] += pow(10, -12);
       }
   }

   int nrows = q[0]->basis->numFuncs;
   int ncols = q[0]->z.len;
//
//   ////////////////////////////////////////////////////////////
//   // serial norm and times
   CMatrix JSerial = CMatrix_init(nrows, ncols);
   printf("\n");
   printf("timing serial GetJacobian\n");
   for(int i = 0; i < numTests; ++i)
   {
      TIME(GetJacobian(q[i], JSerial));
      printf("Frobenius norm of GetJacobian = %.16e\n", CMatrixFrobenius(JSerial));
      printf("\n");
   }
   ////////////////////////////////////////////////////////////


   ////////////////////////////////////////////////////////////
   // parallel norm and times
   CMatrix JParallel = CMatrix_init(nrows, ncols);
   for(int i = 0; i < numTests; ++i)
      QuadAllocBasisOmp(q[i], omp_get_max_threads());

   printf("timing parallel GetJacobianOmp\n");
   for(int i = 0; i < numTests; ++i)
   {
      TIME(GetJacobianOmp(q[i], JParallel));
      printf("Frobenius norm of GetJacobianOmp = %.16e\n", CMatrixFrobenius(JParallel));
      printf("\n");
   }
   printf("\n");
   ////////////////////////////////////////////////////////////

   printf("maximum absolute difference = %.16e\n", CMatrixMaxDifference(JSerial, JParallel));

   CMatrix_free(JSerial);
   CMatrix_free(JParallel);
   for(int i = 0; i < numTests; ++i)
      QuadFreeBasisOmp(q[i]);
   for(int i = 0; i < numTests; ++i)
      quadrature_free(q[i]);

   return EXIT_SUCCESS;
}
