/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */

#include "NodeElimination.h"

#include "GetJacobian.h"
#include "LeastSquaresNewton.h"
#include "ConstrainedOptimization.h"
#include "Quadrature.h"
#include "LINALG.h"
#include "Print.h"
#include "Conditional_Debug.h"
#include "get_time.h"

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <omp.h>
#include "../plasma-17.1/include/plasma.h"

extern int MAX_DIM;
static void FreeMemory(quadrature *q_temp, quadrature *q_new);
RMatrix PredictorLapack(quadrature *q, int *arrayIndex);
RMatrix PredictorMergeLapack(quadrature *q, int *arrayIndex, int *arrayAssociate);
#ifdef _OPENMP
RMatrix PredictorPlasma(quadrature *q, int *arrayIndex);
#endif
static void ReorderWithBoundaryDist(RMatrix predictor, quadrature *q, int *arrayIndex);
static void ReorderWithBoundaryDistMerge(RMatrix predictor, quadrature *q, int *arrayIndex, int *arrayAssociate);
static void ExtractFromPredictor(RMatrix Z, int arrayIndex, quadrature *q);
static void ExtractFromPredictorMerge(RMatrix Z, int arrayIndex, int arrayAssociate, quadrature *q);
//static void ConstrainQTemp(quadrature *q_new, int arrayIndex, quadrature *q_temp);
static void InsertionSort(int num_entries, double *norms, int *arrayIndex);
static double TwoNorm(int n, double *z);
ATTR_UNUSED static bool TestQR(const CMatrix Q);


/***************************************************************************************************
 * A new node elimination scheme that eliminates one node at a time and computes
 * the initial guess for constrained Newton's method. Subsequently, Newton's method is called
 * to obtain quadrature rule with fewer nodes. The procedure is repeated
 * until no more nodes can be eliminated.
 ***************************************************************************************************/
double PREDICTOR_TIME = 0.0;
void NodeElimination(const quadrature *q_initial, quadrature *q_final, history *hist)
{
   if(q_initial->isFullyInitialized != GQ_TRUE) {
      PRINT_ERR(STR_QUAD_NOT_FULL_INIT, __LINE__, __FILE__);
      return;
   }
   if(q_final->isFullyInitialized != GQ_TRUE) {
      PRINT_ERR(STR_QUAD_NOT_FULL_INIT, __LINE__, __FILE__);
      return;
   }
   int nodesInitial = q_initial->num_nodes;
   int numFuncs     = q_initial->basis->numFuncs;
   int dim          = q_initial->dim;
   double tol       = QUAD_TOL; // 10^(-14);

   RMatrix(*Predictor_Ptr)(quadrature *, int *);
#ifdef _OPENMP
   if(PLASMA_CONDITION()) Predictor_Ptr = &PredictorPlasma;
   else                   Predictor_Ptr = &PredictorLapack;
#else
   Predictor_Ptr = &PredictorLapack;
#endif

   quadrature *q_temp = quadrature_make_full_copy(q_initial);
   quadrature *q_new  = quadrature_make_full_copy(q_initial);

   // test accuracy of the initial quadrature
   // if residual is too large, attempt to find quadrature using Newton's method.
   // if Newton's method finds a solution, proceed to Node elimination, otherwise return.
   double res = QuadTestIntegral(q_initial, orthogonal);
   PrintDouble(res, "initial orthogonal basis residual in NodeElimination");
   if(fabs(res) > tol)
   {
      int its = 0;
      bool SOL_FLAG = LeastSquaresNewton(ON, q_temp, &its);
      if(SOL_FLAG == SOL_FOUND) {
         nodesInitial = q_temp->num_nodes;
         quadrature_realloc_array(nodesInitial, q_new);
         quadrature_assign(q_temp, q_new);
      }
      else if(SOL_FLAG == SOL_NOT_FOUND) {
         Print("Initial quadrature did not converge. The initial guess should be more accurate.\n");
         FreeMemory(q_temp, q_new);
         return;
      }
   }

   // run Node Elimination Algorithm. Theoretical optimum is reached when
   // (dim+1)*k = numFuncs, at which point elimination is pursued no further.
   int n_prev            = -1;
   int n_cur             = nodesInitial;
   bool_enum CONSTR_FLAG = OFF;
   bool SOL_FLAG         = SOL_NOT_FOUND;
   int n_opt             = ceil(1.0*numFuncs/(dim+1));
   double efficiency     = (double)n_opt/nodesInitial;
   PrintElimInfo( dim, n_cur , n_opt, efficiency);
   while( ((dim+1)*n_cur > numFuncs)  && (n_cur >= 1) )
   {
      if(dim == MAX_DIM)
         if(n_cur == n_prev) CONSTR_FLAG = ON;
         else                CONSTR_FLAG = OFF;
      else CONSTR_FLAG = OFF;
      n_prev = n_cur;

      quadrature_shrink_array(n_cur-1, q_temp);
      quadrature_shrink_array(n_cur, q_new);

//      int *arrayIndexMerge = (int *)malloc(SIZE_INT(n_cur*(n_cur-1)/2));
//      int *arrayAssociate = (int *)malloc(SIZE_INT(n_cur*(n_cur-1)/2));
//      RMatrix ZMerge = PredictorMergeLapack(q_new, arrayIndexMerge, arrayAssociate);
//      for(int i = 0; i < n_cur*(n_cur-1)/2; ++i)
//      {
//         // store ith initial quadrature guess of Z in q_temp
//         ExtractFromPredictorMerge(ZMerge, arrayIndexMerge[i], arrayAssociate[arrayIndexMerge[i]], q_temp);
//         if(V_InfNorm(q_temp->z) >= QUAD_HUGE)
//            continue;
//
//         int its = 0;
//         SOL_FLAG = LeastSquaresNewton(CONSTR_FLAG, q_temp, &its);
//         // store nodes and weights if Newton's method succeeded, update history
//         if(SOL_FLAG == SOL_FOUND)
//         {
//            n_cur = q_temp->num_nodes;
//            efficiency = (double)n_opt/n_cur;
//            quadrature_realloc_array(q_temp->num_nodes, q_new);
//            quadrature_assign(q_temp, q_new);
//
//            hist->hist_array[hist->total_elims].nodes_tot = n_cur;
//            hist->hist_array[hist->total_elims].success_node = i;
//            hist->hist_array[hist->total_elims].success_its = its;
//            ++hist->total_elims;
//            break; // break i loop
//         }
//         else if(SOL_FLAG == SOL_NOT_FOUND)
//            quadrature_realloc_array(n_cur-1, q_temp); // replace with a cheaper routine
//      }// end i loop
//      RMatrix_free(ZMerge);
//      free(arrayIndexMerge);
//      free(arrayAssociate);

      int *arrayIndex = (int *)malloc(SIZE_INT(n_cur));
      double start_time = get_cur_time();
      RMatrix Z = Predictor_Ptr(q_new, arrayIndex);
      PREDICTOR_TIME += get_cur_time() - start_time;
      ReorderWithBoundaryDist(Z, q_temp, arrayIndex);
      for(int i = 0; i < n_cur; ++i)
      {
         // store ith initial quadrature guess of Z in q_temp
         ExtractFromPredictor(Z, arrayIndex[i], q_temp);
         if(V_InfNorm(q_temp->z) >= QUAD_HUGE)
            continue;

         int its = 0;
         SOL_FLAG = LeastSquaresNewton(CONSTR_FLAG, q_temp, &its);
         // store nodes and weights if Newton's method succeeded, update history
         if(SOL_FLAG == SOL_FOUND)
         {
            n_cur = q_temp->num_nodes;
            efficiency = (double)n_opt/n_cur;
            quadrature_realloc_array(q_temp->num_nodes, q_new);
            quadrature_assign(q_temp, q_new);

            hist->hist_array[hist->total_elims].nodes_tot = n_cur;
            hist->hist_array[hist->total_elims].success_node = i;
            hist->hist_array[hist->total_elims].success_its = its;
            ++hist->total_elims;
            break; // break i loop
         }
         else if(SOL_FLAG == SOL_NOT_FOUND)
            quadrature_realloc_array(n_cur-1, q_temp); // replace with a cheaper routine
      }// end i loop
      RMatrix_free(Z);
      free(arrayIndex);

      if(SOL_FLAG == SOL_FOUND && CONSTR_FLAG == ON)
         Print("Solution found using constraint");
      if (SOL_FLAG == SOL_FOUND)
         PrintElimInfo(dim, n_cur, n_opt, efficiency);
      else if( (SOL_FLAG == SOL_NOT_FOUND) && (CONSTR_FLAG != ON)  && (dim == MAX_DIM) )
         Print("rerunning with constrained optimization");
      // end NodeElimination if repeat flag is off and all iterations were unsuccessful
      else if(SOL_FLAG == SOL_NOT_FOUND) {
         Print("Last iteration did not converge");
         break; // break while loop
      }
   }// end while loop

   // save nodes and weights
   q_final->num_nodes = n_cur;
   quadrature_assign(q_new, q_final);
   res = QuadTestIntegral(q_final, orthogonal);
   PrintDouble(res, "Final orthogonal basis residual in NodeElimination"); printf("\n");

   FreeMemory(q_temp, q_new);
}// end NodeElimination


static void FreeMemory(quadrature *q_temp, quadrature *q_new)
{
   quadrature_free(q_temp);
   quadrature_free(q_new);
}


RMatrix PredictorLapack(quadrature *q, int *arrayIndex)
{
   const int numFuncs  = q->basis->numFuncs;
   const int dim       = q->dim;
   const int n_cur     = q->num_nodes;
   const int nrows     = numFuncs;
   const int ncols     = (dim+1) * n_cur;
   assert(nrows < ncols);

   CMatrix J = CMatrix_init(nrows, ncols);
   GetJacobian(q, J);
   CMatrix J_TR = CMatrix_Transpose(J);

   // construct QR factorization of transpose of the jacobian
   int N_REFL   = MIN(J_TR.rows, J_TR.cols);
   Vector REFL = Vector_init(N_REFL);
   if(DGEQRF_LAPACK(J_TR, REFL) != 0)
      PRINT_ERR(STR_LAPACK_ERR, __LINE__, __FILE__);

   RMatrix Z           = RMatrix_init(n_cur, n_cur*(dim+1));
   RMatrix dZ          = RMatrix_init(n_cur, n_cur*(dim+1));
   CMatrix QFull       = CMatrix_init(n_cur*(dim+1), n_cur*(dim+1));
   CMatrix Q2          = CMatrix_init(n_cur*(dim+1), n_cur*(dim+1)-numFuncs);
   CMatrix Q2Weight    = CMatrix_init(n_cur, n_cur*(dim+1)-numFuncs);
   CMatrix Q2Mult      = CMatrix_init(n_cur*(dim+1), n_cur);

   // obtain Q(from J_TR) explicitly(multiply by identity)
   for(int i = 0; i < QFull.rows; ++i) C_ELEM_ID(QFull, i, i) = 1.0;
   if(DORMQR_LAPACK('L',  'N', REFL, J_TR, QFull) != 0)
      PRINT_ERR(STR_LAPACK_ERR, __LINE__, __FILE__);

   // extract Q2
   for(int i = 0; i < Q2.rows; ++i)
      for(int j = 0; j < Q2.cols; ++j)
         C_ELEM_ID(Q2, i, j) = C_ELEM_ID(QFull, i, j+numFuncs);
   // extract ith rows of Q2 corresponding to weights
   for(int i = 0; i < Q2Weight.rows; ++i)
      for(int j = 0; j < Q2Weight.cols; ++j)
         C_ELEM_ID(Q2Weight, i, j) = C_ELEM_ID(QFull, i, j+numFuncs);

   double *ith_row = (double *)malloc(Q2Weight.cols*sizeof(double));
   double *norm_Q2_ROW = (double *)malloc(n_cur*sizeof(double));
   for(int i = 0; i < n_cur; ++i) {
      for(int j = 0; j < Q2Weight.cols; ++j)
         ith_row[j] = C_ELEM_ID(Q2Weight, i, j);
      norm_Q2_ROW[i] = TwoNorm(Q2Weight.cols, ith_row);
   }

   // multiply Q2 by its weights rows, product goes to columns
   Q2Weight = CMatrix_Transpose(Q2Weight);
   DGEMM_LAPACK(Q2, Q2Weight, Q2Mult);

   for(int i = 0; i < Q2Mult.cols; ++i)
      for(int j = 0; j < Q2Mult.rows; ++j)
         dZ.rid[i][j] = C_ELEM_ID(Q2Mult, j, i) * q->w[i]/SQUARE(norm_Q2_ROW[i]);

   for(int i = 0; i < Z.rows; ++i)
      for(int j = 0; j < Z.cols; ++j)
         Z.rid[i][j] = q->z.id[j] - dZ.rid[i][j];

   double *distance  = (double *)malloc(SIZE_DOUBLE(dZ.rows));
   for(int i = 0; i < dZ.rows; ++i) {
      distance[i] = TwoNorm(ncols, dZ.rid[i]);
      arrayIndex[i] = i;
   }
   InsertionSort(n_cur, distance, arrayIndex);
//   PrintDoubles(distance, n_cur, "distance");

   free(ith_row);
   free(norm_Q2_ROW);
   free(distance);
   Vector_free(REFL);
   CMatrix_free(J_TR);
   RMatrix_free(dZ);
   CMatrix_free(QFull);
   CMatrix_free(Q2);
   CMatrix_free(Q2Weight);
   CMatrix_free(Q2Mult);

   return Z;
}

#ifdef _OPENMP
RMatrix PredictorPlasma(quadrature *q, int *arrayIndex)
{
   const int numFuncs  = q->basis->numFuncs;
   const int dim       = q->dim;
   const int n_cur     = q->num_nodes;
   const int nrows     = numFuncs;
   const int ncols     = (dim+1) * n_cur;
   int INFO;
   assert(nrows < ncols);

   RMatrix Z           = RMatrix_init(n_cur, n_cur*(dim+1));
   RMatrix dZ          = RMatrix_init(n_cur, n_cur*(dim+1));
   CMatrix QFull       = CMatrix_init(n_cur*(dim+1), n_cur*(dim+1));
   CMatrix Q2          = CMatrix_init(n_cur*(dim+1), n_cur*(dim+1)-numFuncs);
   CMatrix Q2Weight    = CMatrix_init(n_cur, n_cur*(dim+1)-numFuncs);
   CMatrix Q2Mult      = CMatrix_init(n_cur*(dim+1), n_cur);
   CMatrix J           = CMatrix_init(numFuncs, (dim+1)*n_cur);

   plasma_init();

   // construct QR factorization of transpose of the jacobian
   GetJacobian(q, J);
   CMatrix J_TR = CMatrix_Transpose(J);
   int LDJ = J_TR.rows;
   plasma_desc_t T;
   INFO = plasma_dgeqrf(J_TR.rows, J_TR.cols, J_TR.id, LDJ, &T);
   if(INFO != 0)
      PRINT_ERR(STR_PLASMA_ERR, __LINE__, __FILE__);

   // obtain Q(from J_TR) explicitly(multiply by identity)
   for(int i = 0; i < QFull.rows; ++i) C_ELEM_ID(QFull, i, i) = 1.0;
   int n_refl = MIN(J_TR.rows, J_TR.cols);
   INFO = plasma_dormqr(PlasmaRight, PlasmaNoTrans,
                 QFull.rows, QFull.cols,
                 n_refl, J_TR.id, LDJ, T,
                 QFull.id, QFull.rows);
   if(INFO != 0)
      PRINT_ERR(STR_PLASMA_ERR, __LINE__, __FILE__);

   // extract Q2
   for(int i = 0; i < Q2.rows; ++i)
      for(int j = 0; j < Q2.cols; ++j)
         C_ELEM_ID(Q2, i, j) = C_ELEM_ID(QFull, i, j+numFuncs);
   // extract ith rows of Q2 corresponding to weights
   for(int i = 0; i < Q2Weight.rows; ++i)
      for(int j = 0; j < Q2Weight.cols; ++j)
         C_ELEM_ID(Q2Weight, i, j) = C_ELEM_ID(QFull, i, j+numFuncs);

   double *ith_row = (double *)malloc(Q2Weight.cols*sizeof(double));
   double *norm_Q2_ROW = (double *)malloc(n_cur*sizeof(double));
   for(int i = 0; i < n_cur; ++i) {
      for(int j = 0; j < Q2Weight.cols; ++j)
         ith_row[j] = C_ELEM_ID(Q2Weight, i, j);
      norm_Q2_ROW[i] = TwoNorm(Q2Weight.cols, ith_row);
   }

   // multiply Q2 by its weights rows, product goes to columns
   Q2Weight = CMatrix_Transpose(Q2Weight);
   DGEMM_PLASMA(Q2, Q2Weight, Q2Mult);

   plasma_desc_destroy(&T);
   plasma_finalize();

   for(int i = 0; i < Q2Mult.cols; ++i)
      for(int j = 0; j < Q2Mult.rows; ++j)
         dZ.rid[i][j] = C_ELEM_ID(Q2Mult, j, i) * q->w[i]/SQUARE(norm_Q2_ROW[i]);

   for(int i = 0; i < Z.rows; ++i)
      for(int j = 0; j < Z.cols; ++j)
         Z.rid[i][j] = q->z.id[j] - dZ.rid[i][j];

   double *distance  = (double *)malloc(SIZE_DOUBLE(dZ.rows));
   for(int i = 0; i < dZ.rows; ++i) {
      distance[i] = TwoNorm(ncols, dZ.rid[i]);
      arrayIndex[i] = i;
   }
   InsertionSort(n_cur, distance, arrayIndex);
//   PrintDoubles(distance, n_cur, "distance");

   free(ith_row);
   free(norm_Q2_ROW);
   free(distance);
   CMatrix_free(J_TR);
   RMatrix_free(dZ);
   CMatrix_free(QFull);
   CMatrix_free(Q2);
   CMatrix_free(Q2Weight);
   CMatrix_free(Q2Mult);

   return Z;
}
#endif

RMatrix PredictorMergeLapack(quadrature *q, int *arrayIndex, int *arrayAssociate)
{
   const int numFuncs  = q->basis->numFuncs;
   const int dim       = q->dim;
   const int n_cur     = q->num_nodes;
   const int nrows     = numFuncs;
   const int ncols     = (dim+1) * n_cur;
   assert(nrows < ncols);

   CMatrix J = CMatrix_init(nrows, ncols);
   GetJacobian(q, J);
   CMatrix J_TR = CMatrix_Transpose(J);
   // construct QR factorization of transpose of the jacobian
   int N_REFL   = MIN(J_TR.rows, J_TR.cols);
   Vector REFL = Vector_init(N_REFL);
   if(DGEQRF_LAPACK(J_TR, REFL) != 0)
      PRINT_ERR(STR_LAPACK_ERR, __LINE__, __FILE__);

   CMatrix MergeMatrix_K  = CMatrix_init(dim, n_cur*(dim+1)-numFuncs);
   CMatrix MergeMatrix_L  = CMatrix_init(dim, n_cur*(dim+1)-numFuncs);
   CMatrix MergeMatrix_KL = CMatrix_init(dim, n_cur*(dim+1)-numFuncs);
   Vector rhs             = Vector_init( MAX(dim, n_cur*(dim+1)-numFuncs) );
   CMatrix Q2             = CMatrix_init(n_cur*(dim+1), n_cur*(dim+1)-numFuncs);
   CMatrix FullQ          = CMatrix_init(n_cur*(dim+1), n_cur*(dim+1));
   RMatrix Z              = RMatrix_init(n_cur*(n_cur-1)/2, n_cur*(dim+1));
   RMatrix dZ             = RMatrix_init(n_cur*(n_cur-1)/2, n_cur*(dim+1));
   double *distance       = (double *)malloc(SIZE_DOUBLE(n_cur*(n_cur-1)/2));
   double *norm_Q2_ROW    = (double *)malloc(n_cur*(n_cur-1)/2*sizeof(double));
   double *dz             = (double *)malloc(n_cur*(dim+1)*sizeof(double));

   //extract full Q and Q2; it will be identical for all k and l
   for(int i = 0; i < FullQ.rows; ++i) C_ELEM_ID(FullQ, i, i) = 1.0;
   if(DORMQR_LAPACK('R',  'N', REFL, J_TR, FullQ) != 0)
      PRINT_ERR(STR_LAPACK_ERR, __LINE__, __FILE__);
   for(int i = 0; i < Q2.rows; ++i)
      for(int j = 0; j < Q2.cols; ++j)
         C_ELEM_ID(Q2, i, j) = C_ELEM_ID(FullQ, i, j+numFuncs);

   int count = 0;
   for(int k = 0; k < n_cur; ++k)
   {
      for(int l = k+1; l < n_cur; ++l)
      {
         // extract info for kth node
         for(int i = 0; i < MergeMatrix_K.rows; ++i)
            for(int j = 0; j < MergeMatrix_K.cols; ++j)
               C_ELEM_ID(MergeMatrix_K, i, j) = C_ELEM_ID(Q2, n_cur+k*dim+i, j);

         // extract info for lth node
         for(int i = 0; i < MergeMatrix_L.rows; ++i)
            for(int j = 0; j < MergeMatrix_L.cols; ++j)
               C_ELEM_ID(MergeMatrix_L, i, j) = C_ELEM_ID(Q2, n_cur+l*dim+i, j);

         // compute merge matrix for q_kl
         for(int i = 0; i < MergeMatrix_KL.rows; ++i)
            for(int j = 0; j < MergeMatrix_KL.cols; ++j)
               C_ELEM_ID(MergeMatrix_KL, i, j) = C_ELEM_ID(MergeMatrix_K, i, j) - C_ELEM_ID(MergeMatrix_L, i, j);

         // compute merge right-hand side for q_kl
         double xkl[MergeMatrix_KL.rows];
         for(int i = 0; i < MergeMatrix_KL.rows; ++i)
            xkl[i] = q->x[dim*l+i] - q->x[dim*k+i];

         //solve dz for q_kl
         memset(rhs.id, 0, rhs.len*sizeof(double));
         for(int i = 0; i < MergeMatrix_KL.rows; ++i)
            rhs.id[i] = xkl[i];
         DGELS_LAPACK(MergeMatrix_KL, rhs);

         // Store new quadrature for q_kl in (kl) row of Z. Rows of Q2 are accessed in the order that
         // corresponds to the order in which quadrature stores its nodes and weights.
         memset(dz, 0, n_cur*(dim+1)*sizeof(double));
         for(int i = 0; i < Q2.rows; ++i) {
            for(int j = 0; j < Q2.cols; ++j)
               dz[i] += C_ELEM_ID(Q2, i , j) * rhs.id[j];
         }

         for(int i = 0; i < q->z.len; ++i)
            R_ELEM_ID(Z, count, i) = q->z.id[i];
         for(int i = 0; i < q->z.len; ++i)
            R_ELEM_ID(Z, count, i) += dz[i];

         R_ELEM_ID(Z, count, k) += R_ELEM_ID(Z, count, l); // add lth weight to k-th weight

         distance[count] = TwoNorm(ncols, dz);
         arrayAssociate[count] = l;
         ++count;
      }
   }

   InsertionSort(n_cur*(n_cur-1)/2, distance, arrayIndex);
//   PrintDoubles(distance, n_cur, "distance");
//   printf("best arrayIndex in merge = %i\n", arrayIndex[0]);
//   printf("node to be eliminated = %i\n", arrayAssociate[arrayIndex[0]]);
//   PrintNodesAndWeights(q, "q_prev");

   Vector_free(REFL);
   CMatrix_free(J_TR);
   CMatrix_free(MergeMatrix_K);
   CMatrix_free(MergeMatrix_L);
   CMatrix_free(MergeMatrix_KL);
   Vector_free(rhs);
   CMatrix_free(Q2);
   CMatrix_free(FullQ);
   RMatrix_free(dZ);
   free(distance);
   free(norm_Q2_ROW);
   free(dz);

   return Z;
}

static void ReorderWithBoundaryDist(RMatrix predictor, quadrature *q, int *arrayIndex)
{
   int numGuesses = predictor.rows;
   quadrature *q_temp = quadrature_make_full_copy(q);
   int maxBound = 5;

   int rowBoundOrder[numGuesses];
   double boundIndex[numGuesses];
   int rowRemainOrder[numGuesses];

   int boundCount = 0;
   int remainCount = 0;
   for(int i = 0; i < numGuesses; ++i)
   {
      ExtractFromPredictor(predictor, arrayIndex[i], q_temp);
      if(QuadInConstraint(q_temp) && boundCount < maxBound) {
         rowBoundOrder[boundCount] = arrayIndex[i];
         boundIndex[boundCount] = QuadMinDistFromTheBoundary(q_temp);
         ++boundCount;
      }
      else {
         rowRemainOrder[remainCount] = arrayIndex[i];
         ++remainCount;
      }
   }

   int tempBoundReorder[numGuesses];
   InsertionSort(boundCount, boundIndex, tempBoundReorder);
   int boundDecreaseReorder[boundCount];
   for(int i = 0; i < boundCount; ++i)
      boundDecreaseReorder[i] = tempBoundReorder[boundCount-1-i];

   for(int i = 0; i < boundCount; ++i)
      arrayIndex[i] = rowBoundOrder[boundDecreaseReorder[i]];
   for(int i = 0; i < remainCount; ++i)
      arrayIndex[boundCount+i] = rowRemainOrder[i];

   quadrature_free(q_temp);
}

static void ReorderWithBoundaryDistMerge(RMatrix predictor, quadrature *q, int *arrayIndex, int *arrayAssociate)
{
   int numGuesses = predictor.rows;
   quadrature *q_temp = quadrature_make_full_copy(q);

   int rowBoundOrder[numGuesses];
   double boundIndex[numGuesses];
   int rowRemainOrder[numGuesses];

   int boundCount = 0;
   int remainCount = 0;
   for(int i = 0; i < numGuesses; ++i)
   {
      ExtractFromPredictorMerge(predictor, arrayIndex[i], arrayAssociate[arrayIndex[i]], q_temp);
      if(QuadInConstraint(q_temp)) {
         rowBoundOrder[boundCount] = arrayIndex[i];
         boundIndex[boundCount] = QuadMinDistFromTheBoundary(q_temp);
         ++boundCount;
      }
      else {
         rowRemainOrder[remainCount] = arrayIndex[i];
         ++remainCount;
      }
   }

   int tempBoundReorder[numGuesses];
   InsertionSort(boundCount, boundIndex, tempBoundReorder);
   int boundDecreaseReorder[boundCount];
   for(int i = 0; i < boundCount; ++i)
      boundDecreaseReorder[i] = tempBoundReorder[boundCount-1-i];

   for(int i = 0; i < boundCount; ++i)
      arrayIndex[i] = rowBoundOrder[boundDecreaseReorder[i]];
   for(int i = 0; i < remainCount; ++i)
      arrayIndex[boundCount+i] = rowRemainOrder[i];

   quadrature_free(q_temp);
}

static void ExtractFromPredictor(RMatrix Z, int arrayIndex, quadrature *q)
{
   int count, j, d;
   int dim = q->dim;

   for(count = 0, j = 0; j < Z.cols/(dim+1); ++j)
   {
      if(j == arrayIndex) continue;

      q->w[count] = Z.rid[arrayIndex][j];
      for(d = 0; d < dim; ++d)
         q->x[count*dim+d] = Z.rid[arrayIndex][Z.cols/(dim+1)+j*dim+d];
      ++count;
   }
}

static void ExtractFromPredictorMerge(RMatrix Z, int arrayIndex, int arrayAssociate, quadrature *q)
{
   int count, j, d;
   int dim = q->dim;

   for(count = 0, j = 0; j < Z.cols/(dim+1); ++j)
   {
      if(j == arrayAssociate) continue;

      q->w[count] = Z.rid[arrayIndex][j];
      for(d = 0; d < dim; ++d)
         q->x[count*dim+d] = Z.rid[arrayIndex][Z.cols/(dim+1)+j*dim+d];
      ++count;
   }
}

static void InsertionSort(int num_entries, double *norms, int *arrayIndex)
{
   for(int i = 0; i < num_entries; ++i)
   {
      int j = i;
      double temp = norms[i];
      while( (j > 0) && (norms[j-1] > temp) ) {
         arrayIndex[j] = arrayIndex[j-1];
         norms[j] = norms[j-1];
         --j;
      }
      arrayIndex[j] = i;
      norms[j] = temp;
   }
}

static double TwoNorm(int n, double *z)
{
   assert(n >= 1);

   double norm = 0.0;
   for(int i = 0; i < n; ++i)
      norm += SQUARE(z[i]);

   norm = SQRT(norm);
   return norm;
}

// assumes Q has less rows than columns
static bool TestQR(const CMatrix Q)
{
   int rows = Q.rows; int cols = Q.cols;
   double max_non_diag = 0.0, max_diag = 1.0;
   double non_diag = 0.0, diag = 0.0;
   double temp = -1.0, tol = POW_DOUBLE(10, -12);

   CMatrix Q_COPY = CMatrix_init(rows, cols);
   memcpy(Q_COPY.id, Q.id, SIZE_DOUBLE(Q.len));

   CMatrix Q_COPY_TR = CMatrix_init(cols, rows);
   Transpose(cols, rows, Q_COPY.id, Q_COPY_TR.id);

   CMatrix Identity = CMatrix_init(rows, rows);
   DGEMM_LAPACK(Q_COPY, Q_COPY_TR, Identity);

   for(int i = 0; i < rows; ++i)
   {
      temp = fabs(Identity.id[i*rows+i]);
      diag = MAX(temp, diag);

      for(int j = 0; j < rows; ++j)
      {
         if(j == i) continue;
         temp = fabs(Identity.id[i*rows+j]);
         non_diag = MAX(temp, non_diag);
      }
   }

   CMatrix_free(Q_COPY);
   CMatrix_free(Q_COPY_TR);
   CMatrix_free(Identity);

   double err1 = fabs(diag - max_diag);
   double err2 = fabs(non_diag - max_non_diag);
   double max_error = MAX(err1, err2);

   if(max_error > tol) return FAILED;
   else                return PASSED;
}


