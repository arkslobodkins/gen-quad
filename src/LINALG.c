#include "LINALG.h"
#include "Print.h"

#include <stdlib.h>
#include <omp.h>
#include <plasma.h>


int DGEMM_LAPACK(CMatrix A, CMatrix B, CMatrix C)
{
   char TRANS1  = 'N';
   char TRANS2  = 'N';
   double alpha = 1.0;
   double beta  = 0.0;
   int M = A.rows;
   int K = A.cols;
   int N = B.cols;
   if(A.rows != C.rows) return INV_INPUT;
   if(A.cols != B.rows) return INV_INPUT;
   if(B.cols != C.cols) return INV_INPUT;

   dgemm_(&TRANS1, &TRANS2, &M, &N, &K, &alpha, A.id, &M, B.id, &K, &beta, C.id, &M);

   return GQ_SUCCESS;
}

int DGEQR2_LAPACK(CMatrix A, Vector TAU)
{
   assert(TAU.len == MIN(A.rows, A.cols));

   int INFO;
   int LDA = A.rows;
   double *WORK = (double *)malloc(SIZE_DOUBLE(2*A.cols));
   dgeqr2_(&A.rows, &A.cols, A.id, &LDA, TAU.id, WORK, &INFO);
   free(WORK);

   return INFO;
}

int DGEQRF_LAPACK(CMatrix A, Vector TAU)
{
   assert(TAU.len == MIN(A.rows, A.cols));

   int INFO;
   int LDA = A.rows;
   int LWORK = 2*A.cols;
   double *WORK = (double *)malloc(SIZE_DOUBLE(LWORK));
   dgeqrf_(&A.rows, &A.cols, A.id, &LDA, TAU.id, WORK, &LWORK, &INFO);
   free(WORK);

   return INFO;
}

int DORMQR_LAPACK(char SIDE, char TRANS, Vector TAU, CMatrix Q, CMatrix A)
{
   int INFO;
   int LDA = A.rows;
   int LDQ = Q.rows;
   int lworkQR = SIDE == 'L' ? 2*A.cols : 2*A.rows;
   double *workQR = (double *)malloc(SIZE_DOUBLE(lworkQR));
   dormqr_(&SIDE, &TRANS, &A.rows, &A.cols, &TAU.len, Q.id, &LDQ, TAU.id, A.id, &LDA, workQR, &lworkQR, &INFO);
   free(workQR);

   return INFO;
}

int DORGQR_LAPACK(CMatrix Q, Vector TAU)
{
   int INFO;
   int LDA = Q.rows;
   int LWORK = 2*Q.cols;
   double *WORK = (double *)malloc(LWORK*sizeof(LWORK));

   dorgqr_(&Q.rows, &Q.cols, &TAU.len, Q.id, &LDA,
          TAU.id, WORK, &LWORK, &INFO);

   free(WORK);
   return INFO;
}

int DGESVD_LAPACK(CMatrix A, Vector VMin)
{
   char JOBU  = 'N'; // U is not computed
   char JOBVT = 'S'; //  return MIN(M, N) rows of V^T, i.e. right singular vectors

   int M = A.rows;
   int N = A.cols;
   int LDA = A.rows;
   Vector SINGV = Vector_init(MIN(M, N));
   if(N > M)
      PRINT_WARN("svd received underdetermined matrix", __LINE__, __FILE__);

   double *U = NULL;
   int LDU = 1;
   int LDVT = MIN(M, N);               // expected to be N in general
   CMatrix VT = CMatrix_init(LDVT, N); // expected to be N x N in general

   int LWORK = 8*MIN(M, N); // assumes M is larger than N, 5*MIN is the minimum work required
   Vector WORK = Vector_init(LWORK);
   int INFO;

   dgesvd_(&JOBU, &JOBVT, &M, &N, A.id, &LDA,
           SINGV.id, U, &LDU, VT.id, &LDVT,
           WORK.id, &LWORK, &INFO);
//   printf("minimum singular value = %.12e\n", SINGV.id[MIN(M, N)-1]);
   CMatrix_GetRow(MIN(M, N)-1, VT, VMin); // get last row of VT;

   Vector_free(SINGV);
   Vector_free(WORK);
   CMatrix_free(VT);

   if(!INFO) return LAPACK_ERR;
   else      return GQ_SUCCESS;
}

int DGELS_LAPACK(CMatrix A, Vector RHS_TO_X)
{
   assert(RHS_TO_X.len == MAX(A.rows, A.cols));

   int INFO;
   char TRANS = 'N';
   int NRHS = 1;
   int LDA = A.rows;
   int LEAD_DIM = MAX(A.rows, A.cols);
   int LWORK = 5*A.cols;
   Vector WORK = Vector_init(LWORK);

   dgels_(&TRANS, &A.rows, &A.cols, &NRHS, A.id, &LDA,
          RHS_TO_X.id, &LEAD_DIM, WORK.id, &WORK.len, &INFO);
   Vector_free(WORK);

   return INFO;
}

#ifdef _OPENMP
int DGELS_PLASMA(CMatrix A, Vector RHS_TO_X)
{
   assert(RHS_TO_X.len == MAX(A.rows, A.cols));

   int NRHS = 1;
   int LDA = A.rows;
   int LEAD_DIM = MAX(A.rows, A.cols);

   plasma_init();
   plasma_desc_t T;
   int INFO = plasma_dgels(PlasmaNoTrans,
                           A.rows, A.cols, NRHS,
                           A.id, LDA, &T,
                           RHS_TO_X.id, LEAD_DIM);
   plasma_desc_destroy(&T);
   plasma_finalize();

   return INFO;
}

int DGEMM_PLASMA(CMatrix A, CMatrix B, CMatrix C)
{
   plasma_init();
   char TRANS1  = PlasmaNoTrans;
   char TRANS2  = PlasmaNoTrans;
   double alpha = 1.0;
   double beta  = 0.0;
   int M = A.rows;
   int K = A.cols;
   int N = B.cols;
   int INFO = plasma_dgemm(TRANS1, TRANS2, M, N, K, alpha, A.id, M,
                           B.id, K, beta, C.id, M);
   plasma_finalize();
   return INFO;
}
#endif

// Simple transpose
void Transpose(int M, int N, const double *A, double *B)
{
   #ifdef _OPENMP
   #pragma omp parallel for default(shared) schedule(static) num_threads(omp_get_max_threads())
   #endif
   for(int i = 0; i < M; ++i)
      for(int j = 0; j < N; ++j)
         B[i+j*M] = A[j+i*N];
}
