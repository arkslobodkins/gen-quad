#include "LINALG.h"
#include "GENERAL_QUADRATURE.h"

#include <stdlib.h>
#include <omp.h>
#include "../plasma-17.1/include/plasma.h"



void MATMUL_LAPACK(int M, int K, int N, double *MAT1, double *MAT2, double *MAT3)
{
   char TRANS1  = 'N';
   char TRANS2  = 'N';
   double alpha = 1.0;
   double beta  = 0.0;
   dgemm_(&TRANS1, &TRANS2, &M, &N, &K, &alpha, MAT1, &M, MAT2, &K, &beta, MAT3, &M);
}


int DGEQR2_LAPACK(CMatrix A, Vector TAU)
{
   assert(TAU.len = MIN(A.rows, A.cols));

   int INFO;
   int LDA = A.rows;
   double *WORK = (double *)malloc(A.cols*size_double);

   dgeqr2_(&A.rows, &A.cols, A.id, &LDA, TAU.id, WORK, &INFO);

   free(WORK);
   return INFO;
}


int DORMQR_LAPACK_M(char SIDE, char TRANS, Vector TAU, CMatrix Q, CMatrix A)
{
   int INFO;
   int LDQ     = SIDE == 'L' ?  Q.rows : Q.cols;
   int lworkQR = SIDE == 'L' ? 2*A.cols : 2*A.rows;
   double *workQR = (double *)malloc(SIZE_DOUBLE(lworkQR));
   dormqr_(&SIDE, &TRANS, &A.rows, &A.cols, &TAU.len, Q.id, &LDQ, TAU.id, A.id, &LDQ, workQR, &lworkQR, &INFO);

   free(workQR);
   return INFO;
}


int DORMQR_LAPACK_V(char SIDE, char TRANS, double *REFL_Q, CMatrix Q, Vector A)
{
   int INFO;
   int LDJ        = Q.rows;
   int ONE_COLUMN = 1;
   int N_REFL     = MIN(Q.rows, Q.cols);
   int lworkQR    = 2*Q.cols;
   double *workQR = (double *)malloc(SIZE_DOUBLE(lworkQR));
   dormqr_(&SIDE, &TRANS, &A.len, &ONE_COLUMN, &N_REFL, Q.id, &LDJ, REFL_Q, A.id, &LDJ, workQR, &lworkQR, &INFO);

   free(workQR);
   return INFO;
}


int DGEQRF_PLASMA(CMatrix A)
{
   int LDA = A.rows;

   plasma_init(omp_get_max_threads());
   plasma_desc_t T;

   int INFO = plasma_dgeqrf(A.rows, A.cols, A.id, LDA, &T);

   plasma_desc_destroy(&T);
   plasma_finalize();

   return INFO;
}


int DORMQR_PLASMA_M(char SIDE, char TRANS, CMatrix Q, CMatrix A)
{
   plasma_enum_t side;
   plasma_enum_t trans;
   int LDA;

   if (SIDE == 'L')     side = PlasmaLeft;
   else if(SIDE == 'R') side = PlasmaRight;
   else                 return -1;

   if(TRANS == 'N')      trans = PlasmaNoTrans;
   else if(TRANS == 'R') trans = PlasmaTrans;
   else                  return -1;

   if (SIDE == 'L')     LDA = A.rows;
   else if(SIDE == 'R') LDA = A.cols;

   int n_refl = MIN(Q.rows, Q.cols);

   plasma_init(omp_get_max_threads());
   plasma_desc_t T;

   int INFO = plasma_dormqr(side, trans,
         A.rows, A.cols,
         n_refl, Q.id, Q.rows,
         T,
         A.id, LDA);

   plasma_desc_destroy(&T);
   plasma_finalize();

   return INFO;
}


int DGELS_LAPACK(CMatrix A, Vector RHS_TO_X)
{
   int LWORK = 5*A.cols;
   Vector WORK = Vector_init(LWORK);

   int INFO;
   char TRANS = 'N';
   int NRHS = 1;
   int LDA = A.rows;
   int LEAD_DIM = MAX(A.rows, A.cols);

   dgels_(&TRANS, &A.rows, &A.cols, &NRHS, A.id, &LDA,
         RHS_TO_X.id, &LEAD_DIM, WORK.id, &WORK.len, &INFO);

   Vector_free(WORK);
   return INFO;
}


int DGELS_PLASMA(CMatrix A, Vector RHS_TO_X)
{
   int NRHS = 1;
   int LDA = A.rows;
   int LEAD_DIM = MAX(A.rows, A.cols);

   plasma_init(omp_get_max_threads());
   plasma_desc_t T;

   int INFO = plasma_dgels(PlasmaNoTrans,
         A.rows, A.cols, NRHS,
         A.id, LDA, &T,
         RHS_TO_X.id, LEAD_DIM);

   plasma_desc_destroy(&T);
   plasma_finalize();

   return INFO;
}

// Simple transpose
void Transpose(int M, int N, const double *A, double *B)
{
   for(int i = 0; i < M; ++i)
      for(int j = 0; j < N; ++j)
         B[i+j*M] = A[j+i*N];
}
