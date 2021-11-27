/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */

#include "NodeElimination.h"

#include "GetJacobian.h"
#include "LeastSquaresNewton.h"
#include "ConstrainedOptimization.h"
#include "Quadrature.h"
#include "Matrix.h"
#include "BasisIndices.h"
#include "LINALG.h"
#include "Print.h"
#include "GENERAL_QUADRATURE.h"
#include "Conditional_Debug.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <assert.h>
#include <stdint.h>
#include <string.h>

extern int MAX_DIM;
static void FreeMemory(INT_8 *basis_indices, quadrature *q_temp, quadrature *q_new);
static RMatrix Predictor(quadrature *q, INT_8 *basis_indices, int *arrayIndex);
static void ExtractFromPredictor(RMatrix Z, int arrayIndex, quadrature *q);
static void ConstrainQTemp(quadrature *q_new, int arrayIndex, quadrature *q_temp);
static void InsertionSort(int num_entries, double *norms, int *arrayIndex);
static double TwoNorm(int n, double *z);
ATTR_UNUSED static bool TestQR(int numRows, int qiCols, const double *Q);


/***************************************************************************************************
 * A new node elimination scheme that eliminates one node at a time and computes
 * the initial guess for constrained Newton's method. Subsequently, Newton's method is called
 * to obtain quadrature rule with fewer nodes. The procedure is repeated
 * until no more nodes can be eliminated.
 ***************************************************************************************************/
void NodeElimination(const quadrature *q_initial, quadrature *q_final, glist *history)
{
   if(q_initial->setFuncsConstrFlag!= 1)
      PRINT_ERR("Functions and constraints for q_initial are not set", __LINE__, __FILE__);
   if(q_final->setFuncsConstrFlag != 1)
      PRINT_ERR("Functions and constraints for q_final are not set",  __LINE__, __FILE__);

   int n_initial       = q_initial->num_nodes;
   int num_funcs       = q_initial->num_funcs;
   int deg             = q_initial->deg;
   int dim             = q_initial->dim;
   int *dims           = q_initial->dims;
   double tol = QUAD_TOL; // 10^(-15);
   LibraryType libType = LAPACK;

   INT_8 *basis_indices = (INT_8 *)malloc( (num_funcs*dim)*sizeof(INT_8) );
   BasisIndices(deg, dim, basis_indices);

   quadrature *q_temp = quadrature_make_full_copy(q_initial);
   quadrature *q_new  = quadrature_make_full_copy(q_initial);

   // test accuracy of the initial quadrature
   // if residual is too large, attempt to find quadrature using Newton's method.
   // if Newton's method finds a solution, proceed to Node elimination, otherwise return.
   double res = QuadTestIntegral(q_initial);
   PrintDouble(res, "initial residual in NodeElimination");
   if(fabs(res) > tol)
   {
      bool SOL_FLAG = LeastSquaresNewton(libType, ON, basis_indices, q_temp, 0);
      if(SOL_FLAG == SOL_FOUND) {
         n_initial = q_temp->num_nodes;
         quadrature_realloc(q_temp->num_nodes, dim, dims, deg, q_new);
         quadrature_assign(q_temp, q_new);
      }
      else if(SOL_FLAG == SOL_NOT_FOUND) {
         Print("Initial quadrature did not converge. The initial guess should be more accurate.\n");
         FreeMemory(basis_indices, q_temp, q_new);
         return;
      }
   }

   // run Node Elimination Algorithm. Theoretical optimum is reached when
   // (dim+1)*k = num_funcs, at which point elimination is pursued no further.
   int n_prev                   = -1;
   int n_cur                    = n_initial;
   bool_enum CONSTR_FLAG        = OFF;
   bool SOL_FLAG                = SOL_NOT_FOUND;
   int n_opt                    = ceil(1.0*num_funcs/(dim+1));
   double efficiency            = (double)n_opt/n_initial;
   PrintElimInfo( dim, n_cur , n_opt, efficiency);
   while( ((dim+1)*n_cur > num_funcs)  && (n_cur >= 1) )
   {
      if(n_cur <= 1) goto FREERETURN;
      if(dim == MAX_DIM)
         if(n_cur == n_prev) CONSTR_FLAG = ON;
         else                CONSTR_FLAG = OFF;
      else CONSTR_FLAG = OFF;
      n_prev = n_cur;

      // Initialize parameters and arrays
      quadrature_reinit(n_cur-1, q_temp);
      quadrature_reinit(n_cur, q_new);

      // extract initial guesses and run Newton's method
      int *arrayIndex = (int *)malloc(SIZE_INT(n_cur));
      RMatrix Z = Predictor(q_new, basis_indices, arrayIndex);
      for(int i = 0; i < n_cur; ++i)
      {
         // store ith initial quadrature guess of Z in q_temp
         ExtractFromPredictor(Z, arrayIndex[i], q_temp);
         if(V_InfNorm(q_temp->z) >= QUAD_HUGE) {
            char error_string[80];
            strcpy(error_string, STR_QUAD_HUGE_ERR);
            strcat(error_string, ", this is expected to happen sometimes");
            PRINT_ERR(error_string, __LINE__, __FILE__);
            continue;
         }

         int its = 0;
         ConstrainQTemp(q_new, arrayIndex[i], q_temp);
         SOL_FLAG = LeastSquaresNewton(libType, CONSTR_FLAG, basis_indices, q_temp, &its);
         // store nodes and weights if Newton's method succeeded, update history
         if(SOL_FLAG == SOL_FOUND)
         {
            n_cur = q_temp->num_nodes;
            efficiency = (double)n_opt/n_cur;
            quadrature_realloc(q_temp->num_nodes, dim, dims, deg, q_new);
            quadrature_assign(q_temp, q_new);

            hist_data *hist_d = (hist_data *)malloc(sizeof(hist_data));
            hist_d->nodes_tot = n_cur;
            hist_d->success_node = i;
            hist_d->success_its = its;
            glist_push(history, (void *)hist_d);
            break; // break i loop
         }
         else if(SOL_FLAG == SOL_NOT_FOUND)
            quadrature_realloc(n_cur-1, dim, dims, deg, q_temp); // replace with a cheaper routine
      }// end i loop
      RMatrix_free(Z);
      free(arrayIndex);

      if(SOL_FLAG == SOL_FOUND && CONSTR_FLAG == ON)
         Print("Solution found using constraint");

      // end NodeElimination if repeat flag is off and all iterations were unsuccessful
      if (SOL_FLAG == SOL_FOUND)
         PrintElimInfo(dim, n_cur, n_opt, efficiency);
      else if( (SOL_FLAG == SOL_NOT_FOUND) && (CONSTR_FLAG != ON)  && (dim == MAX_DIM) )
         Print("rerunning with constrained optimization");
      else if(SOL_FLAG == SOL_NOT_FOUND) {
         Print("Last iteration did not converge");
         break; // break while loop
      }

   }// end while loop

FREERETURN:
   // save nodes and weights
   q_final->num_nodes = n_cur;
   quadrature_assign(q_new, q_final);
   res = QuadTestIntegral(q_final);
   PrintDouble(res, "Final residual in NodeElimination"); printf("\n");

   FreeMemory(basis_indices, q_temp, q_new);
}// end NodeElimination


static void FreeMemory(INT_8 *basis_indices, quadrature *q_temp, quadrature *q_new)
{
   free(basis_indices);
   quadrature_free(q_temp);
   quadrature_free(q_new);
}

// not an easy world we live in
static RMatrix Predictor(quadrature *q, INT_8 *basis_indices, int *arrayIndex)
{
   int num_funcs = q->num_funcs;
   int dim = q->dim;
   int n_cur = q->num_nodes;
   const int nrows = num_funcs;
   const int ncols = (dim+1) * n_cur;
   CMatrix J = CMatrix_init(nrows, ncols);
   GetJacobian(basis_indices, q, J);
   CMatrix J_TR = CMatrix_Transpose(J);

   // construct QR factorization of transpose of the jacobian
   int N_REFL   = MIN(J_TR.rows, J_TR.cols);
   Vector REFL = Vector_init(N_REFL);
   if(DGEQR2_LAPACK(J_TR, REFL) != 0)
      PRINT_ERR(STR_LAPACK_ERR, __LINE__, __FILE__);

   // initialize i*(dim+1)th rows of Q, which correspond to weight rows.
   // Entries that correspond to weight indices are initialized to 1.
   CMatrix QWEIGHT  = CMatrix_init((dim+1)*n_cur, n_cur);
   for(int i = 0; i < QWEIGHT.cols; ++i) C_ELEM_ID(QWEIGHT, i*(dim+1), i) = 1.0;

   // obtain i*(dim+1)th rows of Q(from J_TR) explicitly and store them in QWEIGHT columns
   if(DORMQR_LAPACK_M('L',  'T', REFL, J_TR, QWEIGHT) != 0)
      PRINT_ERR(STR_LAPACK_ERR, __LINE__, __FILE__);
   COND_TEST_4;
   QWEIGHT = CMatrix_Transpose(QWEIGHT);

   RMatrix Z           = RMatrix_init(QWEIGHT.rows, QWEIGHT.cols);
   RMatrix dZ          = RMatrix_init(QWEIGHT.rows, QWEIGHT.cols);
   CMatrix Q2          = CMatrix_init(QWEIGHT.rows, QWEIGHT.cols);
   double *signif_ind  = (double *)malloc(SIZE_DOUBLE(n_cur));
   double *norm_Q2_ROW = (double *)malloc(QWEIGHT.rows*sizeof(double));

   // extract i*(dim+1)th rows of Q2 and compute its norm, where Q = [Q1, Q2]
   for(int i = 0; i < QWEIGHT.rows; ++i) {
      for(int j = 0; j < num_funcs; ++j)            C_ELEM_ID(Q2, i, j) = 0.0;
      for(int j = num_funcs; j < QWEIGHT.cols; ++j) C_ELEM_ID(Q2, i, j) = C_ELEM_ID(QWEIGHT, i, j);
   }
   double *jth_row = (double *)malloc(QWEIGHT.cols*size_double);
   for(int i = 0; i < QWEIGHT.rows; ++i) {
      for(int j = 0; j < QWEIGHT.cols; ++j)
         jth_row[j] = C_ELEM_ID(Q2, i, j);
      norm_Q2_ROW[i] = TwoNorm(QWEIGHT.cols, jth_row);
   }

   // multiply Q by i*(dim+1)th row of Q2
   Q2 = CMatrix_Transpose(Q2);
   if(DORMQR_LAPACK_M('L', 'N', REFL, J_TR, Q2) != 0)
      PRINT_ERR(STR_LAPACK_ERR, __LINE__, __FILE__);
   Q2 = CMatrix_Transpose(Q2);

   // compute initial guesses for Newton's method and store them in Z
   for(int i = 0; i < QWEIGHT.rows; ++i)
   {
      for(int j = 0; j < QWEIGHT.cols; ++j)
         dZ.rid[i][j] = C_ELEM_ID(Q2, i, j) * q->w[i]/SQUARE(norm_Q2_ROW[i]);
      signif_ind[i] = TwoNorm(ncols, dZ.rid[i]);

      // compute new initial guesses and store them in predictor Z
      for(int j = 0; j < Z.cols/(dim+1); ++j) {
         Z.rid[i][j] = q->w[j] - dZ.rid[i][(dim+1)*j];
         for(int d = 0; d < dim; ++d)
            Z.rid[i][n_cur+j*dim+d] = q->x[j*dim+d] - dZ.rid[i][j*(dim+1)+d+1];
      }
   }
   InsertionSort(n_cur, signif_ind, arrayIndex); // store indices of signif_ind in arrayIndex in ascending order

   free(jth_row);
   free(norm_Q2_ROW);
   free(signif_ind);
   Vector_free(REFL);
   RMatrix_free(dZ);
   CMatrix_free(J_TR);
   CMatrix_free(QWEIGHT);
   CMatrix_free(Q2);

   return Z;
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

static void ConstrainQTemp(quadrature *q_new, int arrayIndex, quadrature *q_temp)
{
   ConstrVectData cVectData = ConstrVectDataInit();
   if(QuadInConstraint(q_temp) == false) {
      quadrature *q_prev = quadrature_without_element(q_new, arrayIndex);
      ConstrainedOptimization(q_prev, q_temp, &cVectData);

      if(cVectData.N_OR_W == WEIGHT)
         quadrature_remove_element(cVectData.boundaryNodeId, q_temp);
      quadrature_free(q_prev);
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

static bool TestQR(int rows, int cols, const double *Q)
{

   double max_non_diag = 0.0, max_diag = 1.0;
   double non_diag = 0.0, diag = 0.0;
   double temp = -1.0, tol = POW_DOUBLE(10, -12);

   double *Q_COPY = (double *)calloc( rows*cols, sizeof(double) );
   memcpy(Q_COPY, Q, SIZE_DOUBLE(rows*cols));

   double *Q_COPY_TR= (double *)calloc( cols*rows, sizeof(double) );
   Transpose(cols, rows, Q_COPY, Q_COPY_TR);

   double *I = (double *)calloc( cols*cols, sizeof(double) );
   MATMUL_LAPACK(cols, rows, cols, Q_COPY_TR, Q_COPY, I);

   for(int i = 0; i < cols; ++i)
   {
      temp = fabs(I[i*cols+i]);
      diag = MAX(temp, diag);

      for(int j = 0; j < cols; ++j)
      {
         if(j == i) continue;
         temp = fabs(I[i*cols+j]);
         non_diag = MAX(temp, non_diag);
      }
   }

   free(I);
   free(Q_COPY);
   free(Q_COPY_TR);

   double err1 = fabs(diag - max_diag);
   double err2 = fabs(non_diag - max_non_diag);
   double max_error = MAX(err1, err2);

   if(max_error > tol) return FAILED;
   else                return PASSED;
}

