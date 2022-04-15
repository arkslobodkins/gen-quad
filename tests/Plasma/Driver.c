#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <mkl_lapacke.h>
#include <mkl.h>
#include <plasma.h>

#include "get_time.h"

double compute_norm(int len, double *v)
{
   double norm = 0.0;
   for(int i = 0; i < len; ++i)
      norm += v[i]*v[i];

   return sqrt(norm);
}

int main(int argc, char *argv[])
{
   int nrows = 5000;
   int ncols = 9000;
   int NRHS  = 1;
   int LDA   = nrows;
   int LEAD_DIM = nrows>ncols ? nrows:ncols;

   double *A_plasma = (double *)malloc(nrows*ncols*sizeof(double));
   double *b_plasma = (double *)malloc(LEAD_DIM*sizeof(double));

   for(int i = 0; i < nrows*ncols; ++i)
         A_plasma[i] = (double)rand() / (double) RAND_MAX;
   for(int i = 0; i < LEAD_DIM; ++i)
      b_plasma[i] = (double)rand() / (double) RAND_MAX;

   // Allocate copies so that identical problem will be solved with mkl lapack
   double *A_lapack = (double *)malloc(nrows*ncols*sizeof(double));
   memcpy(A_lapack, A_plasma, nrows*ncols*sizeof(double));
   double *b_lapack = (double *)malloc(LEAD_DIM*sizeof(double));
   memcpy(b_lapack, b_plasma, LEAD_DIM*sizeof(double));

   mkl_set_threading_layer(MKL_THREADING_SEQUENTIAL);
   printf("OPENMP enanbled with %i threads\n\n", omp_get_max_threads());

   //////////////////////////////////////////////////////////////////////
   plasma_init(omp_get_max_threads());
   plasma_desc_t T;

   double PL_START = get_cur_time();
   plasma_dgels(PlasmaNoTrans,
                nrows, ncols, NRHS,
                A_plasma, LDA, &T,
                b_plasma, LEAD_DIM);
   double PL_END = get_cur_time();
   double pl_sol_norm = compute_norm(ncols, b_plasma);

   free(A_plasma);
   free(b_plasma);
   plasma_desc_destroy(&T);
   plasma_finalize();

   printf("plasma time = %lf\n", PL_END-PL_START);
   printf("plasma sol norm = %.16e\n", pl_sol_norm);
   //////////////////////////////////////////////////////////////////////


   //////////////////////////////////////////////////////////////////////
   double LP_START = get_cur_time();
   LAPACKE_dgels(LAPACK_COL_MAJOR, 'N', nrows, ncols, NRHS,
                 A_lapack, LDA, b_lapack, LEAD_DIM);
   double LP_END = get_cur_time();
   double la_sol_norm = compute_norm(ncols, b_lapack);

   free(A_lapack);
   free(b_lapack);

   printf("mkl lapack time = %lf\n", LP_END-LP_START);
   printf("mkl lapack sol norm = %.16e\n", la_sol_norm);
   //////////////////////////////////////////////////////////////////////

   return EXIT_SUCCESS;
}
