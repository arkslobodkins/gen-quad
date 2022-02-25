/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */

#include "NodeElimination.h"

#include "GetJacobian.h"
#include "GetFunction.h"
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
#include <plasma.h>

extern int MAX_DIM;

#define SEARCH_DIM 3
#define MAX_FAILS_LEVEL_1 4
#define MAX_FAILS_LEVEL_2 6
#define MAX_FAILS_ELIM 12

#define PASSED true
#define FAILED false

typedef struct
{
   int index;
   double d;
} DistanceStruct;

typedef struct
{
   double t[SEARCH_DIM];
} ShortenParams;
static const ShortenParams shortParams = {{0.5, 0.25, 0.125}};

static DistanceStruct* distance_init(int n)
{
   DistanceStruct *distanceStr = (DistanceStruct *)malloc(n*sizeof(DistanceStruct));
   for(int i = 0; i < n; ++i)  distanceStr[i].index = i;
   for(int i = 0; i < n; ++i)  distanceStr[i].d = 0.0;
   return distanceStr;
}

static bool LsqSearch(bool_enum CONSTR_FLAG, quadrature *q_new, history *hist);
static bool TreeSearch(bool_enum CONSTR_FLAG, quadrature *q_new, history *hist);

static RMatrix PredictorLapack(quadrature *q, DistanceStruct *distance);
#ifdef _OPENMP
static RMatrix PredictorPlasma(quadrature *q, DistanceStruct *distance);
#endif
__attribute__unused static RMatrix PredictorSimple(quadrature *q, DistanceStruct *signifIndex);

static void ExtractFromPredictor(RMatrix Z, int arrayIndex, quadrature *q);
static void ExtractFromPredictorFull(RMatrix Z, int arrayIndex, quadrature *q);
__attribute__unused static void ReorderWithBoundaryDist(RMatrix predictor, quadrature *q, DistanceStruct *arrayIndex);
__attribute__unused static void InsertionSort(int num_entries, double *norms, int *arrayIndex);

static int compareDouble(const void *a, const void *b);
static double TwoNorm(int n, double *z);

__attribute__unused static bool TestQR(const CMatrix Q);
__attribute__unused static void QuadSavePlots(quadrature *q);


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
   quadrature *q_new  = quadrature_make_full_copy(q_initial);

   // test accuracy of the initial quadrature
   // if residual is too large, attempt to find quadrature using Newton's method.
   // if Newton's method finds a solution, proceed to Node elimination, otherwise return.
   double res = QuadTestIntegral(q_initial, orthogonal);
   PrintDouble(res, "initial orthogonal basis residual in NodeElimination");
   if(fabs(res) > tol)
   {
      quadrature *q_temp = quadrature_make_full_copy(q_initial);
      LSQ_out lsq_out = LeastSquaresNewton(ON, q_temp);
      if(lsq_out.SOL_FLAG == SOL_FOUND) {
         quadrature_assign_resize(q_temp, q_new);
         quadrature_free(q_temp);
      }
      else if(lsq_out.SOL_FLAG == SOL_NOT_FOUND) {
         Print("Initial quadrature did not converge. The initial guess should be more accurate.\n");
         quadrature_free(q_temp);
         quadrature_free(q_new);
         return;
      }
   }

   // run Node Elimination Algorithm. Theoretical optimum is reached when
   // (dim+1)*k = numFuncs, at which point elimination is pursued no further.
   int n_cur             = nodesInitial;
   bool SOL_FLAG         = SOL_NOT_FOUND;
   int n_opt             = ceil(1.0*numFuncs/(dim+1));
   double efficiency     = (double)n_opt/nodesInitial;
   PrintElimInfo( dim, n_cur , n_opt, efficiency);
   int repeat = 0;
   while( (n_cur > n_opt)  && (n_cur >= 2) )
   {
      SOL_FLAG = LsqSearch(OFF, q_new, hist);        // perform regular search first
      if(SOL_FLAG == SOL_NOT_FOUND)                  // perform deeper search
      {
         if(repeat == 0)                             // if failed early
            SOL_FLAG = TreeSearch(ON, q_new, hist);
         else if(dim != MAX_DIM)
            SOL_FLAG = TreeSearch(OFF, q_new, hist); // without boundary constraints
         else if(dim == MAX_DIM)
            SOL_FLAG = TreeSearch(ON, q_new, hist);  // add boundary constraints for the last dimension
      }

      n_cur = q_new->num_nodes;
      efficiency = (double)n_opt/n_cur;
      if (SOL_FLAG == SOL_FOUND)
         PrintElimInfo(dim, n_cur, n_opt, efficiency);
      else if(SOL_FLAG == SOL_NOT_FOUND)
      {
         Print("Last iteration did not converge");
         break; // break while loop
      }
      ++repeat;
   }// end while loop

   // save nodes and weights
   quadrature_assign_resize(q_new, q_final);
   res = QuadTestIntegral(q_final, orthogonal);
   PrintDouble(res, "Final orthogonal basis residual in NodeElimination"); printf("\n");

   quadrature_free(q_new);

}// end NodeElimination


static bool LsqSearch(bool_enum CONSTR_FLAG, quadrature *q_new, history *hist)
{
   RMatrix(*Predictor_Ptr)(quadrature *, DistanceStruct *);
   #ifdef _OPENMP
      if(PLASMA_CONDITION()
         && OMP_CONDITION(q_new->deg, q_new->dim))
            Predictor_Ptr = &PredictorPlasma;
      else  Predictor_Ptr = &PredictorLapack;
   #else
      Predictor_Ptr = &PredictorLapack;
   #endif

   int n_cur = q_new->num_nodes;

   DistanceStruct *distanceWeight = distance_init(n_cur);
   RMatrix Z = Predictor_Ptr(q_new, distanceWeight);

   int failCount = 0;
   LSQ_out lsq_out = {0, SOL_NOT_FOUND};
   quadrature *q_temp = quadrature_make_full_copy(q_new);
   ReorderWithBoundaryDist(Z, q_new, distanceWeight);
   for(int i = 0; i < n_cur; ++i)
   {
      quadrature_realloc_array(n_cur-1, q_temp);
      // store ith initial quadrature guess in q_temp
      ExtractFromPredictor(Z, distanceWeight[i].index, q_temp);

      if(V_InfNorm(q_temp->z) >= QUAD_HUGE) continue;
      if(!QuadInConstraint(q_temp))         continue;

      lsq_out = LeastSquaresNewton(CONSTR_FLAG, q_temp);
      // store nodes and weights if Newton's method succeeded, update history
      if(lsq_out.SOL_FLAG == SOL_FOUND)
      {
         n_cur = q_temp->num_nodes;
         quadrature_assign_resize(q_temp, q_new);

         hist->hist_array[hist->total_elims].nodes_tot    = n_cur;
         hist->hist_array[hist->total_elims].success_node = failCount;
         hist->hist_array[hist->total_elims].success_its  = lsq_out.its;
         ++hist->total_elims;
         break;
      }
      else if(lsq_out.SOL_FLAG == SOL_NOT_FOUND)
      {
         quadrature_realloc_array(n_cur-1, q_temp);
         if(++failCount >= MAX_FAILS_ELIM)
            break;
      }
   }
   quadrature_free(q_temp);
   free(distanceWeight);
   RMatrix_free(Z);

   return lsq_out.SOL_FLAG;
}


static bool TreeSearch(bool_enum CONSTR_FLAG, quadrature *q_new, history *hist)
{
   RMatrix(*Predictor_Ptr)(quadrature *, DistanceStruct *);
   #ifdef _OPENMP
      if(PLASMA_CONDITION()
         && OMP_CONDITION(q_new->deg, q_new->dim))
           Predictor_Ptr = &PredictorPlasma;
      else Predictor_Ptr = &PredictorLapack;
   #else
      Predictor_Ptr = &PredictorLapack;
   #endif

   int n_cur = q_new->num_nodes;
   int dim = q_new->dim;
   bool SOL_FLAG = SOL_NOT_FOUND;

   DistanceStruct *distance__1 = distance_init(n_cur);
   RMatrix Z__1 = Predictor_Ptr(q_new, distance__1);

   quadrature *qnewtemp__1 = quadrature_make_full_copy(q_new);
   quadrature *qsearch__1  = quadrature_make_full_copy(q_new);
   Vector dz__1 = Vector_init(n_cur*(dim+1));
   int failCount__1 = 0;
   for(int i1 = 0; i1 < n_cur; ++i1)
   {

      // level 1 search
      quadrature_assign(q_new, qnewtemp__1);
      ExtractFromPredictorFull(Z__1, distance__1[i1].index, qsearch__1);
      bool searchElimFlag__1 = SOL_NOT_FOUND;
      bool earlyConstraint   = false;
      LSQ_out searchNewGuessFlag__1;
      int sd1;
      for(sd1 = 0; sd1 < SEARCH_DIM; ++sd1)
      {
         VectorAddScale(1.0, qsearch__1->z, -1.0, qnewtemp__1->z, dz__1);
         VectorAddScale(1.0, qnewtemp__1->z, shortParams.t[sd1], dz__1, qsearch__1->z);
         searchNewGuessFlag__1 = LeastSquaresNewton(CONSTR_FLAG, qsearch__1);

         if(qsearch__1->num_nodes != n_cur)
         {
            quadrature_assign_resize(qsearch__1, q_new);
            SOL_FLAG = SOL_FOUND;
            earlyConstraint = true;
            printf("succeeded at %i th iteration, with constraints at depth level 0.5, and damping parameter %.4e\n", failCount__1, shortParams.t[sd1]);
            hist->hist_array[hist->total_elims].nodes_tot    = q_new->num_nodes;
            hist->hist_array[hist->total_elims].success_node = failCount__1;
            hist->hist_array[hist->total_elims].success_its  = searchNewGuessFlag__1.its;
            ++hist->total_elims;
            break; // break sd1 loop
         }


         if(searchNewGuessFlag__1.SOL_FLAG)
            break; // break sd1 loop
      }

      if(earlyConstraint)
         break; // break i1 loop

      if(searchNewGuessFlag__1.SOL_FLAG)
      {
         searchElimFlag__1 = LsqSearch(CONSTR_FLAG, qsearch__1, hist);

         // if eliminated a node, save solution and break all loops
         if(searchElimFlag__1)
         {
            quadrature_assign_resize(qsearch__1, q_new);
            SOL_FLAG = SOL_FOUND;
            printf("succeeded at %i th iteration, with shortening at depth level 1, and damping parameter %.4e\n", failCount__1, shortParams.t[sd1]);
            break; // break i1 loop if new solution was found
         }
         // new initial guess was found, but node was not eliminated, take new initial guess to level 2
         else if(!searchElimFlag__1)
            quadrature_assign(qsearch__1, qnewtemp__1);
      }
      else if(!searchNewGuessFlag__1.SOL_FLAG)
      {
         if(++failCount__1 >= MAX_FAILS_LEVEL_1)
         {
            SOL_FLAG = SOL_NOT_FOUND;
            break;
         }
         continue; // go to next node
      }


      // level 2 search
      // assumes searchNewGuessFlag__1 is true and searchElimFlag__1 false;
      assert(searchNewGuessFlag__1.SOL_FLAG  && !searchElimFlag__1);

      quadrature *qsearch__2 = quadrature_make_full_copy(qnewtemp__1);
      DistanceStruct *distance__2 = distance_init(n_cur);
      Vector dz__2 = Vector_init(n_cur*(dim+1));
      RMatrix Z__2  = Predictor_Ptr(qnewtemp__1, distance__2);

      bool searchElimFlag__2 = SOL_NOT_FOUND;
      int failCount__2 = 0;
      for(int i2 = 0; i2 < n_cur; ++i2)
      {
         ExtractFromPredictorFull(Z__2, distance__2[i2].index, qsearch__2);
         LSQ_out searchNewGuessFlag__2;

         int sd2 = 0;
         for(sd2 = 0; sd2 < SEARCH_DIM; ++sd2)
         {
            VectorAddScale(1.0, qsearch__2->z, -1.0, qnewtemp__1->z, dz__2);
            VectorAddScale(1.0, qnewtemp__1->z, shortParams.t[sd2], dz__2, qsearch__2->z);
            searchNewGuessFlag__2 = LeastSquaresNewton(CONSTR_FLAG, qsearch__2);

            if(qsearch__2->num_nodes != n_cur)
            {
               quadrature_assign_resize(qsearch__2, q_new);
               SOL_FLAG = SOL_FOUND;
               earlyConstraint = true;
               printf("succeeded at %i->%i th iteration, with constraints at depth level 1.5, and damping parameter %.4e\n", failCount__1, failCount__2, shortParams.t[sd1]);
               hist->hist_array[hist->total_elims].nodes_tot    = q_new->num_nodes;
               hist->hist_array[hist->total_elims].success_node = failCount__2;
               hist->hist_array[hist->total_elims].success_its  = searchNewGuessFlag__2.its;
               ++hist->total_elims;
               break; // break sd2 loop
            }

            if(searchNewGuessFlag__2.SOL_FLAG)
               break; // break sd2 loop
         }
         if(earlyConstraint)
            break; // break i2 loop

         if(searchNewGuessFlag__2.SOL_FLAG)
         {
            searchElimFlag__2 = LsqSearch(CONSTR_FLAG, qsearch__2, hist);
            // if eliminated a node, save solution and break all loops
            if(searchElimFlag__2)
            {
               quadrature_assign_resize(qsearch__2, q_new);
               SOL_FLAG = SOL_FOUND;
               printf("succeeded at %i->%i th iteration, with shortening at depth level 2, and damping parameter %.4e\n", failCount__1, failCount__2, shortParams.t[sd2]);
               break; // break i2 loop and later i1 loop
            }
         }

         if(++failCount__2 >= MAX_FAILS_LEVEL_2)
            break;
      }
      quadrature_free(qsearch__2);
      RMatrix_free(Z__2);
      Vector_free(dz__2);
      free(distance__2);

      if(earlyConstraint)
         break;

      if(searchElimFlag__2)
         break;

      if(++failCount__1 >= MAX_FAILS_LEVEL_1)
      {
         SOL_FLAG = SOL_NOT_FOUND;
         break;
      }

   } // end i loop

   quadrature_free(qnewtemp__1);
   quadrature_free(qsearch__1);
   RMatrix_free(Z__1);
   Vector_free(dz__1);
   free(distance__1);
   return SOL_FLAG;
}


// assumes distance indices are initialized properly
static RMatrix PredictorLapack(quadrature *q, DistanceStruct *distance)
{
   const int numFuncs  = q->basis->numFuncs;
   const int dim       = q->dim;
   const int n_cur     = q->num_nodes;
   const int nrows     = numFuncs;
   const int ncols     = (dim+1) * n_cur;
   assert(nrows < ncols);

   CMatrix J = CMatrix_init(nrows, ncols);
   GetJacobian(q, J);
   double start_time__1 = get_cur_time();
   CMatrix J_TR = CMatrix_Transpose(J);

   // construct QR factorization of transpose of the Jacobian
   int N_REFL   = MIN(J_TR.rows, J_TR.cols);
   Vector REFL = Vector_init(N_REFL);
  if(DGEQRF_LAPACK(J_TR, REFL) != 0)
      PRINT_ERR(STR_LAPACK_ERR, __LINE__, __FILE__);

   RMatrix Z        = RMatrix_init(n_cur, n_cur*(dim+1));                  // M x N
   RMatrix dZ       = RMatrix_init(n_cur, n_cur*(dim+1));                  // M x N
   CMatrix QFull    = CMatrix_init(n_cur*(dim+1), n_cur*(dim+1));          // N x N
   CMatrix Q2       = CMatrix_init(n_cur*(dim+1), n_cur*(dim+1)-numFuncs); // N x (N-M)
   CMatrix Q2Weight = CMatrix_init(n_cur, n_cur*(dim+1)-numFuncs);         // M x (N-M)
   CMatrix Q2Mult   = CMatrix_init(n_cur*(dim+1), n_cur);                  // N x M

   // obtain Q(from J_TR) explicitly(multiply by identity)
   for(int i = 0; i < QFull.rows; ++i) C_ELEM_ID(QFull, i, i) = 1.0;
   if(DORMQR_LAPACK('L',  'N', REFL, J_TR, QFull) != 0)
      PRINT_ERR(STR_LAPACK_ERR, __LINE__, __FILE__);

//   COND_TEST_3;

   // loop over columns for more efficient access(although marginal compared to QR)
   // extract Q2
   for(int j = 0; j < Q2.cols; ++j)
      for(int i = 0; i < Q2.rows; ++i)
         C_ELEM_ID(Q2, i, j) = C_ELEM_ID(QFull, i, j+numFuncs);
   // extract ith rows of Q2 corresponding to weights
   for(int j = 0; j < Q2Weight.cols; ++j)
      for(int i = 0; i < Q2Weight.rows; ++i)
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

   int LA_FLAG;
   if((LA_FLAG = DGEMM_LAPACK(Q2, Q2Weight, Q2Mult)) != 0)
      PRINT_ERR(ERR_STRING(LA_FLAG), __LINE__, __FILE__);

   double *w = q->w;
   for(int i = 0; i < Q2Mult.cols; ++i)
      for(int j = 0; j < Q2Mult.rows; ++j)
         dZ.rid[i][j] = C_ELEM_ID(Q2Mult, j, i) * w[i]/SQUARE(norm_Q2_ROW[i]);

   Vector qz = q->z;
   for(int i = 0; i < Z.rows; ++i)
      for(int j = 0; j < Z.cols; ++j)
         Z.rid[i][j] = qz.id[j] - dZ.rid[i][j];

   for(int i = 0; i < dZ.rows; ++i)
      distance[i].d = TwoNorm(ncols, dZ.rid[i]);

   qsort(distance, n_cur, sizeof(DistanceStruct), compareDouble);

   free(ith_row);
   free(norm_Q2_ROW);
   Vector_free(REFL);
   CMatrix_free(J_TR);
   RMatrix_free(dZ);
   CMatrix_free(QFull);
   CMatrix_free(Q2);
   CMatrix_free(Q2Weight);
   CMatrix_free(Q2Mult);

   PREDICTOR_TIME += get_cur_time() - start_time__1;
   return Z;
}


#ifdef _OPENMP
static RMatrix PredictorPlasma(quadrature *q, DistanceStruct *distance)
{
   const int numFuncs  = q->basis->numFuncs;
   const int dim       = q->dim;
   const int n_cur     = q->num_nodes;
   const int nrows     = numFuncs;
   const int ncols     = (dim+1) * n_cur;
   int INFO;
   assert(nrows < ncols);

   RMatrix Z        = RMatrix_init(n_cur, n_cur*(dim+1));
   RMatrix dZ       = RMatrix_init(n_cur, n_cur*(dim+1));
   CMatrix QFull    = CMatrix_init(n_cur*(dim+1), n_cur*(dim+1));
   CMatrix Q2       = CMatrix_init(n_cur*(dim+1), n_cur*(dim+1)-numFuncs);
   CMatrix Q2Weight = CMatrix_init(n_cur, n_cur*(dim+1)-numFuncs);
   CMatrix Q2Mult   = CMatrix_init(n_cur*(dim+1), n_cur);
   CMatrix J        = CMatrix_init(numFuncs, (dim+1)*n_cur);

   // construct QR factorization of transpose of the Jacobian
   QuadAllocBasisOmp(q, omp_get_max_threads());
   GetJacobianOmp(q, J);
   QuadFreeBasisOmp(q, omp_get_max_threads());
   double start_time__1 = get_cur_time();
   CMatrix J_TR = CMatrix_Transpose(J);
   int LDJ = J_TR.rows;
   plasma_init();
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
//   COND_TEST_3;

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
   plasma_desc_destroy(&T);
   plasma_finalize();

   // multiply Q2 by its weights rows, product goes to columns
   Q2Weight = CMatrix_Transpose(Q2Weight);
   DGEMM_PLASMA(Q2, Q2Weight, Q2Mult);

   double *w = q->w;
   for(int i = 0; i < Q2Mult.cols; ++i)
      for(int j = 0; j < Q2Mult.rows; ++j)
         dZ.rid[i][j] = C_ELEM_ID(Q2Mult, j, i) * w[i]/SQUARE(norm_Q2_ROW[i]);

   Vector qz = q->z;
   for(int i = 0; i < Z.rows; ++i)
      for(int j = 0; j < Z.cols; ++j)
         Z.rid[i][j] = qz.id[j] - dZ.rid[i][j];

   for(int i = 0; i < dZ.rows; ++i)
      distance[i].d = TwoNorm(ncols, dZ.rid[i]);

   qsort(distance, n_cur, sizeof(DistanceStruct), compareDouble);

   free(ith_row);
   free(norm_Q2_ROW);
   CMatrix_free(J_TR);
   RMatrix_free(dZ);
   CMatrix_free(QFull);
   CMatrix_free(Q2);
   CMatrix_free(Q2Weight);
   CMatrix_free(Q2Mult);

   PREDICTOR_TIME += get_cur_time() - start_time__1;
   return Z;
}
#endif


static RMatrix PredictorSimple(quadrature *q, DistanceStruct *distance)
{
   const int dim = q->dim;
   const int n_cur = q->num_nodes;
   RMatrix Z = RMatrix_init(n_cur, n_cur*(dim+1));

   CubeParams params;
   params.deg = q->deg;
   params.dim = q->dim;
   CubeParams *paramsPtr = &params;
   BasisInterface interface = SetCubeBasisInterface();
   Basis *cubeBasis = BasisInit((void *)paramsPtr, &interface);

   for(int i = 0; i < Z.rows; ++i)
      for(int j = 0; j < Z.cols; ++j)
         Z.rid[i][j] = q->z.id[j];

   for(int i = 0; i < Z.rows; ++i)
      distance[i].d = 0.0;

   double *x;
   Vector basisFuncs = cubeBasis->functions;
   for(int i = 0; i < Z.rows; ++i)
   {
      x = &q->x[i*dim];
      BasisMonomial(cubeBasis, x, basisFuncs);
      for(int j = 0; j < basisFuncs.len; ++j)
         distance[i].d += SQUARE(basisFuncs.id[j]);
   }
   for(int i = 0; i < Z.rows; ++i)
      distance[i].d *= q->w[i];
   qsort(distance, n_cur, sizeof(DistanceStruct), compareDouble);

   BasisFree(cubeBasis);
   return Z;
}


static void ReorderWithBoundaryDist(RMatrix predictor, quadrature *q, DistanceStruct *distance)
{
   int numGuesses = predictor.rows;
   quadrature *q_temp = quadrature_make_full_copy(q);
   int maxBound = 10;

   double boundDistance[maxBound];
   int boundRowIndex[maxBound];
   int remainRowIndex[numGuesses];
   memset(boundDistance, -1, maxBound*sizeof(double));
   memset(boundRowIndex, -1, maxBound*sizeof(int));
   memset(remainRowIndex, -1, numGuesses*sizeof(int));

   double *tempDistances = (double *)malloc(SIZE_DOUBLE(numGuesses));
   for(int i = 0; i < numGuesses; ++i)
      tempDistances[i] = distance[i].d;

   int boundCount = 0;
   int remainCount = 0;
   for(int i = 0; i < numGuesses; ++i)
   {
      ExtractFromPredictor(predictor, distance[i].index, q_temp);
      if(QuadInConstraint(q_temp) && boundCount < maxBound) {
         boundRowIndex[boundCount] = distance[i].index;
         boundDistance[boundCount] = QuadMinDistFromTheBoundary(q_temp);
         ++boundCount;
      }
      else {
         remainRowIndex[remainCount] = distance[i].index;
         ++remainCount;
      }
   }

   int tempBoundReorder[numGuesses];
   for(int i = 0; i < numGuesses; ++i)
      tempBoundReorder[i] = i;

   InsertionSort(boundCount, boundDistance, tempBoundReorder);

   int boundDecreaseReorder[boundCount];
   double boundDistanceReorder[boundCount];
   for(int i = 0; i < boundCount; ++i) {
      boundDecreaseReorder[i] = tempBoundReorder[boundCount-1-i];
      boundDistanceReorder[i] = boundDistance[boundCount-1-i];
   }

   for(int i = 0; i < boundCount; ++i) {
      distance[i].index = boundRowIndex[boundDecreaseReorder[i]];
      distance[i].d = boundDistanceReorder[i];
   }
   for(int i = 0; i < remainCount; ++i) {
      distance[boundCount+i].index = remainRowIndex[i];
      distance[boundCount+i].d = tempDistances[boundCount+i];
   }
   quadrature_free(q_temp);
   free(tempDistances);
}


static void InsertionSort(int num_entries, double *norms, int *arrayIndex)
{
   for(int i = 0; i < num_entries; ++i)
   {
      int j = i;
      double temp = norms[i];
      while( (j > 0) && (norms[j-1] > temp) )
      {
         arrayIndex[j] = arrayIndex[j-1];
         norms[j] = norms[j-1];
         --j;
      }
      arrayIndex[j] = i;
      norms[j] = temp;
   }
}


// extracts quadrature that belongs to arrayIndex row, exluding the 0 weight
static void ExtractFromPredictor(RMatrix Z, int arrayIndex, quadrature *q)
{
   int count, j, d;
   int dim = q->dim;
   int numNodes = Z.cols/(dim+1);

   for(count = 0, j = 0; j < numNodes; ++j) {
      if(j == arrayIndex) continue;
      q->w[count] = Z.rid[arrayIndex][j];
      ++count;
   }

   for(count = 0, j = 0; j < numNodes; ++j) {
      if(j == arrayIndex) continue;
      for(d = 0; d < dim; ++d)
         q->x[count*dim+d] = Z.rid[arrayIndex][numNodes+j*dim+d];
      ++count;
   }
}


// extracts full quadrature that belongs to arrayIndex row, including the node that contains 0 weight
static void ExtractFromPredictorFull(RMatrix Z, int arrayIndex, quadrature *q)
{
   int dim = q->dim;
   int numNodes = Z.cols/(dim+1);

   for(int j = 0; j < numNodes; ++j)
      q->w[j] = Z.rid[arrayIndex][j];

   for(int j = 0; j < numNodes*dim; ++j)
      q->x[j] = Z.rid[arrayIndex][numNodes+j];
}


// comparison routine to be passed to qsort
static int compareDouble(const void *a, const void *b)
{
   DistanceStruct *ad = (DistanceStruct *)a;
   DistanceStruct *bd = (DistanceStruct *)b;
   if(ad->d > bd->d)
      return 1;
   else return -1;
}


static double TwoNorm(int n, double *z)
{
   assert(n >= 1);
   int incz = 1;
   double norm = dnrm2_(&n, z, &incz);
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


// currently saves evolution of quadrature nodes on a 2-dimensional CUBE(square)
static void QuadSavePlots(quadrature *q)
{
   int dim = q->dim;
   int deg = q->deg;
   int n = q->num_nodes;
   static int plotNum = 0;

   if(dim == 2)
   {
      FILE *gnuplot = popen("gnuplot", "w");
      fprintf(gnuplot, "dim = %i;", (int)dim);
      fprintf(gnuplot, "deg = %i;", (int)deg);
      fprintf(gnuplot, "plotNum = %i;", (int)plotNum);
      fprintf(gnuplot, "numPoints = %i;", (int)n);

      int num_commands = 8;
      char dataPath[80];
      sprintf(dataPath, "../data/quadData%i.txt", plotNum);
      FILE *data = fopen(dataPath, "w");
      for (int i = 0; i < q->num_nodes; i++)
         fprintf(data, "%g %g\n", q->x[2*i], q->x[2*i+1]);

      const char * commandsForGnuplot[] =
      {
         "filePath(n) = sprintf('../data/quadData%i.txt', n)",
         "fileName(a, b, c) = sprintf('../data/data_cube_deg%i_dim%i_#%i.png',a, b, c)",
         "set xrange[0:1]",
         "set yrange [0:1]",
         "set terminal png size 500, 500",
         "set output fileName(deg, dim, plotNum)",
         "set title sprintf('Quadrature on a SQUARE, degree %i, %i points', deg, numPoints)",
         "plot filePath(plotNum) with points pointtype 13 notitle"
      };

      for(int i = 0; i < num_commands; ++i)
         fprintf(gnuplot, "%s \n", commandsForGnuplot[i]);

      fflush(gnuplot);
      fclose(data);
   }
   ++plotNum;
}
