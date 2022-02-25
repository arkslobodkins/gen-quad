#include <stdlib.h>
#include <lapacke.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#include "get_time.h"
#include "../../plasma-17.1/include/plasma.h"

void dgels_(char *TRANS, int *M, int *N, int* nrhs,
            double* A, int* LDA,
            double* B, int* LDB,
            double *WORK, int *LWORK, int *INFO);

int plasma_dgels(plasma_enum_t trans,
                 int m, int n, int nrhs,
                 double *pA, int lda,
                 plasma_desc_t *T,
                 double *pB, int ldb);

double compute_norm(int len, double *v)
{
   double norm = 0.0;
   for(int i = 0; i < len; ++i)
      norm += v[i]*v[i];

   return sqrt(norm);
}

int main(int argc, char *argv[])
{
#ifdef _OPENMP
   printf("OPENMP enanbled with %i threads\n\n", omp_get_max_threads());
#endif

   int nrows = 4000;
   int ncols = 2000;
   int NRHS  = 1;
   int LDA   = nrows;
   int LEAD_DIM = nrows>ncols ? nrows:ncols;


   plasma_init(omp_get_max_threads());

   int retval;
   int IONE = 1;
   int ISEED[4] = {0,0,0,1};
   plasma_desc_t T;

   double *A = (double *)malloc(nrows*ncols*sizeof(double));
   retval = LAPACKE_dlarnv(IONE, ISEED, nrows*ncols, A);
   assert(retval == 0);

   double *b = (double *)malloc(nrows*sizeof(double));
   retval = LAPACKE_dlarnv(IONE, ISEED, LEAD_DIM, b);
   assert(retval == 0);

   for(int i = 0; i < nrows; ++i)
      for(int j = 0; j < ncols; ++j)
         A[j*nrows+i] = (double)rand() / (double) RAND_MAX;

   // Allocate copies so that identical problem will be solved with lapack
   double *A_copy = (double *)malloc(nrows*ncols*sizeof(double));
   memcpy(A_copy, A, nrows*ncols*sizeof(double));

   for(int i = 0; i < nrows; ++i)
      b[i] = (double)rand() / (double) RAND_MAX;

   double *b_copy = (double *)malloc(nrows*sizeof(double));
   memcpy(b_copy, b, nrows*sizeof(double));


   double PL_START = get_cur_time();
   plasma_dgels(PlasmaNoTrans,
                nrows, ncols, NRHS,
                A, LDA,
                &T,
                b, LEAD_DIM);
   double PL_END = get_cur_time();
   double pl_sol_norm = compute_norm(ncols, b);

   free(A);
   free(b);
   plasma_desc_destroy(&T);
   plasma_finalize();

   printf("plasma time = %lf\n", PL_END-PL_START);
   printf("plasma sol norm = %lf\n", pl_sol_norm);




   char TRANS = 'N';
   int INFO;
   int LWORK = 5*nrows;
   double *WORK = (double *)malloc(LWORK*sizeof(double));

   double LP_START = get_cur_time();
   dgels_(&TRANS, &nrows, &ncols, &NRHS,
          A_copy, &LDA,
          b_copy, &LEAD_DIM,
          WORK, &LWORK, &INFO);
   double LP_END = get_cur_time();
   double la_sol_norm = compute_norm(ncols, b_copy);

   free(WORK);
   free(A_copy);
   free(b_copy);

   printf("lapack time = %lf\n", LP_END-LP_START);
   printf("lapack sol norm = %lf\n", la_sol_norm);

   return EXIT_SUCCESS;
}
