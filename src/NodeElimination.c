/* Arkadijs Slobodkins
 * SMU Mathematics
 * Copyright 2022
 */


#include "NodeElimination.h"

#include "GetJacobian.h"
#include "GetFunction.h"
#include "NonlinearSolve.h"
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
#include <mkl_blas.h>

extern int MAX_DIM;

#define SEARCH_DIM 3
#define MAX_FAILS_LEVEL_1 8
#define MAX_FAILS_LEVEL_2 5
#define MAX_FAILS_ELIM 15

#define QR_PASSED true
#define QR_FAILED false

typedef enum { lsq, lev_mar } solver_type;

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


static bool LsqSearch(bool_enum CONSTR_FLAG, quadrature *q_new, history *hist, solver_type st);
static bool WideLsqSearch(bool_enum CONSTR_FLAG, quadrature *q_new, history *hist, solver_type st);
static bool TreeSearch(bool_enum CONSTR_FLAG, quadrature *q_new, history *hist, solver_type st);

static RMatrix PredictorLapack(quadrature *q, DistanceStruct *distance);
#ifdef _OPENMP
static RMatrix PredictorPlasma(quadrature *q, DistanceStruct *distance);
#endif
__attribute__unused static RMatrix PredictorSimple(quadrature *q, DistanceStruct *signifIndex);

static void ExtractFromPredictor(RMatrix Z, int arrayIndex, quadrature *q);
static void ExtractFromPredictorFull(RMatrix Z, int arrayIndex, quadrature *q);
__attribute__unused static int PredictorInConstraintCount(RMatrix Z, const quadrature *q);
__attribute__unused static void ReorderWithBoundaryDist(RMatrix predictor, quadrature *q, DistanceStruct *arrayIndex);
__attribute__unused static void InsertionSort(int num_entries, double *norms, int *arrayIndex);
__attribute__unused static void TestInsertionSort(int num_entries, double *norms); // must be called right after InsertionSort

static int compareDouble(const void *a, const void *b); // comparison routine to be passed to qsort
static double TwoNorm(int n, double *z);

__attribute__unused static bool TestQR(const CMatrix Q);
__attribute__unused static void QuadSavePlots(quadrature *q);


typedef LSQ_out(*Nonlinear_Solve_Ptr)(const bool_enum, quadrature *);
static Nonlinear_Solve_Ptr
set_nonlinear_solve_ptr(solver_type st)
{
   switch(st)
   {
      case lsq:
         return &LeastSquaresNewton;
         break;
      case lev_mar:
         return &LevenbergMarquardt;
      default:
         PRINT_WARN(STR_INV_INPUT, __LINE__, __FILE__);
         return &LeastSquaresNewton;
   }
}


typedef RMatrix(*Predictor_Ptr)(quadrature *, DistanceStruct *);
Predictor_Ptr set_predictor_ptr(quadrature *q)
{
   #ifdef _OPENMP
      if(OMP_CONDITION(q->deg, q->dim))  return &PredictorPlasma;
      else                               return &PredictorLapack;
   #else
      return &PredictorLapack;
   #endif
}

static void UpdateHistory(int n, int snode, int its, int nsols, history *hist)
{
   hist->hist_array[hist->total_elims].nodes_tot    = n;
   hist->hist_array[hist->total_elims].success_node = snode;
   hist->hist_array[hist->total_elims].success_its  = its;
   hist->hist_array[hist->total_elims].num_solutions  = nsols;
   ++hist->total_elims;
}


//////////////////////////////////////////////////////////////////////////////////////////
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


   const int numFuncs = q_initial->basis->numFuncs;
   const int dim      = q_initial->dim;
   const double tol   = QUAD_TOL; // 10^(-14);
   quadrature *q_new  = quadrature_make_full_copy(q_initial);


   // test accuracy of the initial quadrature
   // if residual is too large, attempt to find quadrature using Newton's method.
   // if Newton's method finds a solution, proceed to Node elimination, otherwise return.
   double res = QuadTestIntegral(q_initial, orthogonal);
   PrintDouble(res, "initial orthogonal basis residual in NodeElimination");
   if(fabs(res) > tol)
   {
      quadrature *q_temp = quadrature_make_full_copy(q_initial);
      LSQ_out lsq_out = LevenbergMarquardt(ON, q_temp);
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
   int n_cur             = q_new->num_nodes;
   bool SOL_FLAG         = SOL_NOT_FOUND;
   int n_opt             = ceil(1.0*numFuncs/(dim+1));
   double efficiency     = (double)n_opt/q_new->num_nodes;
   PrintElimInfo(dim, n_cur , n_opt, efficiency);
   while( (n_cur > n_opt) && (n_cur >= 2) )
   {
      SOL_FLAG = LsqSearch(OFF, q_new, hist, lev_mar);
      if(dim == MAX_DIM)
         if(SOL_FLAG == SOL_NOT_FOUND)
         {
            {
               printf("Rerunning with constrained optimization\n");
               SOL_FLAG = LsqSearch(ON, q_new, hist, lev_mar);
            }
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
   }// end while loop
   hist->efficiency = efficiency;

   // save nodes and weights
   quadrature_assign_resize(q_new, q_final);
   res = QuadTestIntegral(q_final, orthogonal);
   PrintDouble(res, "Final orthogonal basis residual in NodeElimination"); printf("\n");

   quadrature_free(q_new);

}// end NodeElimination


static bool LsqSearch(bool_enum CONSTR_FLAG, quadrature *q_new, history *hist, solver_type st)
{
   Predictor_Ptr predictor_ptr             = set_predictor_ptr(q_new);
   Nonlinear_Solve_Ptr nonlinear_solve_ptr = set_nonlinear_solve_ptr(st);

   int n_cur = q_new->num_nodes;
   DistanceStruct *distanceWeight = distance_init(n_cur);
   RMatrix Z = predictor_ptr(q_new, distanceWeight);

   int failCount = 0;
   LSQ_out lsq_out = {0, SOL_NOT_FOUND};
   quadrature *q_temp = quadrature_make_full_copy(q_new);
   quadrature_realloc_array(n_cur-1, q_temp);
   ReorderWithBoundaryDist(Z, q_new, distanceWeight);

   for(int i = 0; i < n_cur; ++i)
   {
      ExtractFromPredictor(Z, distanceWeight[i].index, q_temp); // store ith initial quadrature guess in q_temp

      if(V_InfNorm(q_temp->z) >= QUAD_HUGE) continue;
      if(!QuadInConstraint(q_temp))         continue;

      // update quadrature if Newton's method succeeded, update history
      lsq_out = nonlinear_solve_ptr(CONSTR_FLAG, q_temp);
      if(lsq_out.SOL_FLAG == SOL_FOUND)
      {
         quadrature_assign_resize(q_temp, q_new);
         UpdateHistory(q_temp->num_nodes, failCount, lsq_out.its, 0, hist);
         break;
      }
      else if(lsq_out.SOL_FLAG == SOL_NOT_FOUND)
         if(++failCount >= MAX_FAILS_ELIM)
            break;
   }

   quadrature_free(q_temp);
   free(distanceWeight);
   RMatrix_free(Z);

   return lsq_out.SOL_FLAG;
}




static bool WideLsqSearch(bool_enum CONSTR_FLAG, quadrature *q_new, history *hist, solver_type st)
{
   Predictor_Ptr predictor_ptr             = set_predictor_ptr(q_new);
   Nonlinear_Solve_Ptr nonlinear_solve_ptr = set_nonlinear_solve_ptr(st);

   int n_cur = q_new->num_nodes;

   DistanceStruct *distanceWeight = distance_init(n_cur);
   RMatrix Z = predictor_ptr(q_new, distanceWeight);

   int search_size = MIN(10, q_new->num_nodes);

   int failCount = 0;
   int count = 0;
   LSQ_out lsq_out[search_size];
   quadrature *q_temp[search_size];
   for(int i = 0; i < search_size; ++i)
     q_temp[i] = quadrature_make_full_copy(q_new);

   for(int i = 0; i < n_cur; ++i)
   {
      quadrature_realloc_array(n_cur-1, q_temp[count]);
      ExtractFromPredictor(Z, distanceWeight[i].index, q_temp[count]);

      if(V_InfNorm(q_temp[count]->z) >= QUAD_HUGE) continue;
      if(!QuadInConstraint(q_temp[count]))         continue;

      lsq_out[count] = nonlinear_solve_ptr(CONSTR_FLAG, q_temp[count]);
      if(lsq_out[count].SOL_FLAG == SOL_FOUND)
         ++count;
      else if(lsq_out[count].SOL_FLAG == SOL_NOT_FOUND)
         if(++failCount >= MAX_FAILS_ELIM)
            break;

      if(count >= search_size) break;
   }

   if(count == 0)
   {
      for(int i = 0; i < search_size; ++i)
         quadrature_free(q_temp[i]);
      free(distanceWeight);
      RMatrix_free(Z);
      return SOL_NOT_FOUND;
   }

   Vector min_singv = MinSingvJacobians(count, q_temp);
   VMax max_min_singv = VectorMax(min_singv);
   int successQuad = max_min_singv.max_index;

   n_cur = q_temp[successQuad]->num_nodes;
   quadrature_assign_resize(q_temp[successQuad], q_new);

   UpdateHistory(n_cur, successQuad, lsq_out[successQuad].its, count, hist);

   for(int i = 0; i < search_size; ++i)
      quadrature_free(q_temp[i]);
   free(distanceWeight);
   RMatrix_free(Z);
   Vector_free(min_singv);

   return SOL_FOUND;
}


static bool TreeSearch(bool_enum CONSTR_FLAG, quadrature *q_new, history *hist, solver_type st)
{
   Predictor_Ptr predictor_ptr = set_predictor_ptr(q_new);
   Nonlinear_Solve_Ptr nonlinear_solve_ptr = set_nonlinear_solve_ptr(st);

   int n_cur = q_new->num_nodes;
   int dim = q_new->dim;
   bool SOL_FLAG = SOL_NOT_FOUND;

   DistanceStruct *distance__1 = distance_init(n_cur);
   RMatrix Z__1 = predictor_ptr(q_new, distance__1);

   quadrature *qnewtemp__1 = quadrature_make_full_copy(q_new);
   quadrature *qsearch__1  = quadrature_make_full_copy(q_new);
   quadrature *qstart__1  = quadrature_make_full_copy(q_new);
   Vector dz__1 = Vector_init(n_cur*(dim+1));
   int failCount__1 = 0;
   for(int i1 = 0; i1 < n_cur; ++i1)
   {
      // level 1 search
      ExtractFromPredictorFull(Z__1, distance__1[i1].index, qstart__1);
      if(V_InfNorm(qstart__1->z) >= QUAD_HUGE) continue;
      quadrature_assign(q_new, qnewtemp__1);

      bool searchElimFlag__1 = SOL_NOT_FOUND;
      bool earlyConstraint   = false;
      LSQ_out searchNewGuessFlag__1;
      int sd1;
      // level1 stage1
      for(sd1 = 0; sd1 < SEARCH_DIM; ++sd1)
      {
         VAddScale(1.0, qstart__1->z, -1.0, qnewtemp__1->z, dz__1);                     // dz1 = qstart1 - qnew
         VAddScale(1.0, qnewtemp__1->z, shortParams.t[sd1], dz__1, qsearch__1->z);      // qsearch1 = qnew + ??*dz1
         searchNewGuessFlag__1 = nonlinear_solve_ptr(CONSTR_FLAG, qsearch__1);          // update qsearch1 if success

         // update q_new and break if eliminated node/nodes
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

         // typical case: no nodes eliminated but converged
         if(searchNewGuessFlag__1.SOL_FLAG) break; // break sd1 loop
      }
      // break i1 loop if eliminated node/nodes
      if(earlyConstraint) break;

      // level1 stage 2
      if(searchNewGuessFlag__1.SOL_FLAG)
      {
         assert(qsearch__1->num_nodes == n_cur);
         searchElimFlag__1 = LsqSearch(CONSTR_FLAG, qsearch__1, hist, st);

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
      assert(searchNewGuessFlag__1.SOL_FLAG && !searchElimFlag__1);

      quadrature *qsearch__2 = quadrature_make_full_copy(qnewtemp__1);
      quadrature *qstart__2  = quadrature_make_full_copy(q_new);
      Vector dz__2 = Vector_init(n_cur*(dim+1));
      DistanceStruct *distance__2 = distance_init(n_cur);
      RMatrix Z__2  = predictor_ptr(qnewtemp__1, distance__2);

      bool searchElimFlag__2 = SOL_NOT_FOUND;
      int failCount__2 = 0;
      for(int i2 = 0; i2 < n_cur; ++i2)
      {
         ExtractFromPredictorFull(Z__2, distance__2[i2].index, qstart__2);
         LSQ_out searchNewGuessFlag__2;

         int sd2 = 0;
         // level2 stage1
         for(sd2 = 0; sd2 < SEARCH_DIM; ++sd2)
         {
            VAddScale(1.0, qstart__2->z, -1.0, qnewtemp__1->z, dz__2);                        // dz2 = qstart12 - qnewtemp
            VAddScale(1.0, qnewtemp__1->z, shortParams.t[sd2], dz__2, qsearch__2->z);         // qsearch2 = qnewtemp + ??*dz2
            searchNewGuessFlag__2 = nonlinear_solve_ptr(CONSTR_FLAG, qsearch__2);             // update qsearch2 if success

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


         // level2 stage2
         if(searchNewGuessFlag__2.SOL_FLAG)
         {
            assert(qsearch__2->num_nodes == n_cur);
            searchElimFlag__2 = LsqSearch(CONSTR_FLAG, qsearch__2, hist, st);
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
      quadrature_free(qstart__2);
      RMatrix_free(Z__2);
      Vector_free(dz__2);
      free(distance__2);

      if(earlyConstraint)   break;
      if(searchElimFlag__2) break;

      if(++failCount__1 >= MAX_FAILS_LEVEL_1)
      {
         SOL_FLAG = SOL_NOT_FOUND;
         break;
      }

   } // end i loop

   quadrature_free(qnewtemp__1);
   quadrature_free(qsearch__1);
   quadrature_free(qstart__1);
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
   double start_time__1 = get_cur_time();             // time from here since Jacobian is timed separately
   CMatrix J_TR = CMatrix_Transpose(J, move);

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
   CMatrix_Identity(QFull);
   if(DORMQR_LAPACK('L',  'N', REFL, J_TR, QFull) != 0)
      PRINT_ERR(STR_LAPACK_ERR, __LINE__, __FILE__);

   // extract Q2
   for(int j = 0; j < Q2.cols; ++j)
      #pragma omp simd
      for(int i = 0; i < Q2.rows; ++i)
         C_ELEM_ID(Q2, i, j) = C_ELEM_ID(QFull, i, j+numFuncs);
   // extract ith rows of Q2 corresponding to weights
   for(int j = 0; j < Q2Weight.cols; ++j)
      #pragma omp simd
      for(int i = 0; i < Q2Weight.rows; ++i)
         C_ELEM_ID(Q2Weight, i, j) = C_ELEM_ID(QFull, i, j+numFuncs);

   Vector ith_row = Vector_init(Q2Weight.cols);
   Vector norm_Q2_ROW = Vector_init(Q2Weight.rows);
   for(int i = 0; i < Q2Weight.rows; ++i) {
      CMatrix_GetRow(i, Q2Weight, ith_row);
      norm_Q2_ROW.id[i] = V_TwoNorm(ith_row);
   }

   // multiply Q2 by its weight rows(hence tranpose Q2Weight), product goes to columns
   Q2Weight = CMatrix_Transpose(Q2Weight, move);
   DGEMM_LAPACK(Q2, Q2Weight, Q2Mult);

   double *w = q->w;
   for(int i = 0; i < Q2Mult.cols; ++i)
      #pragma omp simd
      for(int j = 0; j < Q2Mult.rows; ++j)
         dZ.rid[i][j] = C_ELEM_ID(Q2Mult, j, i) * w[i]/SQUARE(norm_Q2_ROW.id[i]);

   Vector qz = q->z;
   for(int i = 0; i < Z.rows; ++i)
      #pragma omp simd
      for(int j = 0; j < Z.cols; ++j)
         Z.rid[i][j] = qz.id[j] - dZ.rid[i][j];

   for(int i = 0; i < dZ.rows; ++i)
      distance[i].d = TwoNorm(dZ.cols, dZ.rid[i]);

   qsort(distance, n_cur, sizeof(DistanceStruct), compareDouble);

   Vector_free(norm_Q2_ROW);
   Vector_free(ith_row);
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
   QuadFreeBasisOmp(q);
   double start_time__1 = get_cur_time();         // time from here since Jacobian is timed separately
   CMatrix J_TR = CMatrix_Transpose(J, move);
   int LDJ = J_TR.rows;
   int n_refl = MIN(J_TR.rows, J_TR.cols);

   plasma_init();
   plasma_desc_t T;
   int INFO = plasma_dgeqrf(J_TR.rows, J_TR.cols, J_TR.id, LDJ, &T);
   if(INFO != 0)
      PRINT_ERR(STR_PLASMA_ERR, __LINE__, __FILE__);

   // obtain Q(from J_TR) explicitly(multiply by identity)
   CMatrix_Identity(QFull);
   INFO = plasma_dormqr(PlasmaRight, PlasmaNoTrans,
                        QFull.rows, QFull.cols,
                        n_refl, J_TR.id, LDJ, T,
                        QFull.id, QFull.rows);
   if(INFO != 0)
      PRINT_ERR(STR_PLASMA_ERR, __LINE__, __FILE__);
   plasma_desc_destroy(&T);
   plasma_finalize();

   // extract Q2
   for(int j = 0; j < Q2.cols; ++j)
      #pragma omp simd
      for(int i = 0; i < Q2.rows; ++i)
         C_ELEM_ID(Q2, i, j) = C_ELEM_ID(QFull, i, j+numFuncs);
   // extract ith rows of Q2 corresponding to weights
   for(int j = 0; j < Q2Weight.cols; ++j)
      #pragma omp simd
      for(int i = 0; i < Q2Weight.rows; ++i)
         C_ELEM_ID(Q2Weight, i, j) = C_ELEM_ID(QFull, i, j+numFuncs);

   Vector ith_row = Vector_init(Q2Weight.cols);
   Vector norm_Q2_ROW = Vector_init(Q2Weight.rows);
   for(int i = 0; i < Q2Weight.rows; ++i) {
      CMatrix_GetRow(i, Q2Weight, ith_row);
      norm_Q2_ROW.id[i] = V_TwoNorm(ith_row);
   }

   // multiply Q2 by its weights rows, product goes to columns
   Q2Weight = CMatrix_Transpose(Q2Weight, move);
   INFO = DGEMM_PLASMA(Q2, Q2Weight, Q2Mult);
   if(INFO != GQ_SUCCESS)
      PRINT_ERR(ERR_STRING(INFO), __LINE__, __FILE__);

   double *w = q->w;
   for(int i = 0; i < Q2Mult.cols; ++i)
      #pragma omp simd
      for(int j = 0; j < Q2Mult.rows; ++j)
         dZ.rid[i][j] = C_ELEM_ID(Q2Mult, j, i) * w[i]/SQUARE(norm_Q2_ROW.id[i]);

   Vector qz = q->z;
   for(int i = 0; i < Z.rows; ++i)
      #pragma omp simd
      for(int j = 0; j < Z.cols; ++j)
         Z.rid[i][j] = qz.id[j] - dZ.rid[i][j];

   for(int i = 0; i < dZ.rows; ++i)
      distance[i].d = TwoNorm(dZ.cols, dZ.rid[i]);

   qsort(distance, n_cur, sizeof(DistanceStruct), compareDouble);

   Vector_free(norm_Q2_ROW);
   Vector_free(ith_row);
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
   double start_time = get_cur_time();

   const int dim = q->dim;
   const int n_cur = q->num_nodes;
   RMatrix Z = RMatrix_init(n_cur, n_cur*(dim+1));

   for(int i = 0; i < Z.rows; ++i)
      for(int j = 0; j < Z.cols; ++j)
         Z.rid[i][j] = q->z.id[j];
   for(int i = 0; i < Z.rows; ++i)
      Z.rid[i][i] = 0.0;

   for(int i = 0; i < Z.rows; ++i)
      distance[i].d = 0.0;

   double *x;
   Vector basisFuncs = q->basis->functions;
   for(int i = 0; i < Z.rows; ++i)
   {
      x = &q->x[i*dim];
      BasisMonomial(q->basis, x, basisFuncs);
      for(int j = 0; j < basisFuncs.len; ++j)
         distance[i].d += SQUARE(basisFuncs.id[j]);
   }
   for(int i = 0; i < Z.rows; ++i)
      distance[i].d *= q->w[i];
   qsort(distance, n_cur, sizeof(DistanceStruct), compareDouble);

   PREDICTOR_TIME += get_cur_time() - start_time;
   return Z;
}

static int PredictorInConstraintCount(RMatrix Z, const quadrature *q)
{
   quadrature *q_temp = quadrature_make_full_copy(q);

   int count = 0;
   for(int i = 0; i < Z.rows; ++i)
   {
      ExtractFromPredictor(Z, i, q_temp);
      if(QuadInConstraint(q_temp))
         ++count;
   }
   quadrature_free(q_temp);

   return count;
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

static void TestInsertionSort(int num_entries, double *norms)
{
   for(int i = 1; i < num_entries; ++i)
      if(norms[i] < norms[i-1])
         PRINT_ERR("InsertionSort test Failed", __LINE__, __FILE__);
}


// extracts quadrature that belongs to arrayIndex row, exluding the 0 weight
static void ExtractFromPredictor(RMatrix Z, int arrayIndex, quadrature *q)
{
   assert(arrayIndex > -1 && arrayIndex < Z.rows);

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
      #pragma omp simd
      for(d = 0; d < dim; ++d)
         q->x[count*dim+d] = Z.rid[arrayIndex][numNodes+j*dim+d];
      ++count;
   }
#ifdef QUAD_DEBUG_ON
   if( fabs(Z.rid[arrayIndex][arrayIndex]) > QUAD_TOL )
      PRINT_ERR("weight should be closer to 0", __LINE__, __FILE__);
#endif
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


static int compareDouble(const void *a, const void *b)
{
   DistanceStruct *ad = (DistanceStruct *)a;
   DistanceStruct *bd = (DistanceStruct *)b;
   return (ad->d > bd->d) - (ad->d < bd->d);
}


static double TwoNorm(int n, double *z)
{
   assert(n >= 1);
   int incz = 1;
   double norm = dnrm2(&n, z, &incz);
   return norm;
}


// assumes Q has greater or equal number of columns than rows
static bool TestQR(const CMatrix Q)
{
   int rows = Q.rows; int cols = Q.cols;
   double max_non_diag = 0.0, max_diag = 1.0;
   double non_diag = 0.0, diag = 0.0;
   double temp = -1.0, tol = POW_DOUBLE(10, -12);

   CMatrix Q_COPY = CMatrix_init(rows, cols);

   CMatrix_Assign(Q, Q_COPY);
   CMatrix Q_COPY_TR = CMatrix_Transpose(Q_COPY, copy);
   CMatrix Identity = CMatrix_init(rows, rows);
   DGEMM_LAPACK(Q_COPY, Q_COPY_TR, Identity);

   for(int i = 0; i < rows; ++i)
   {
      temp = fabs(Identity.id[i*rows+i]);
      diag = MAX(temp, diag);
   }
   for(int i = 0; i < rows; ++i)
      for(int j = 0; j < rows; ++j)
      {
         if(j == i) continue;
         temp = fabs(Identity.id[i*rows+j]);
         non_diag = MAX(temp, non_diag);
      }

   CMatrix_free(Q_COPY);
   CMatrix_free(Q_COPY_TR);
   CMatrix_free(Identity);

   double err1 = fabs(diag - max_diag);
   double err2 = fabs(non_diag - max_non_diag);
   double max_error = MAX(err1, err2);

   if(max_error > tol) return QR_FAILED;
   else                return QR_PASSED;
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
