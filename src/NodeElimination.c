/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */

#include "NodeElimination.h"

#include "GetJacobian.h"
#include "InsertionSort.h"
#include "LeastSquaresNewton.h"
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
static void FreeMemory(int_fast8_t *basis, quadrature *q_temp, quadrature *q_new);
static double TwoNorm(int n, double *z);
static bool TestQR(int numRows, int qiCols, const double *Q);


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

   int i,j,d;
   int n_initial       = q_initial->k;
   int num_funcs       = q_initial->num_funcs;
   int deg             = q_initial->deg;
   int dim             = q_initial->dim;
   int dims[q_initial->num_dims];
   const ATTR_UNUSED DOMAIN_TYPE D = q_initial->D;
   for(int d = 0; d < q_initial->num_dims; ++d)
      dims[d] = q_initial->dims[d];
   double tol = QUAD_TOL; // 10^(-15);

   int_fast8_t *basis = (int_fast8_t *)malloc( (num_funcs*dim)*sizeof(int_fast8_t) );
   BasisIndices(deg, dim, basis);

   quadrature *q_temp = quadrature_make_full_copy(q_initial);
   quadrature *q_new  = quadrature_make_full_copy(q_initial);

   // test accuracy of the initial quadrature
   // if residual is too large, attempt to find quadrature using Newton's method.
   // if Newton's method finds a solution, proceed to Node elimination, otherwise return.
   double res = QuadTestIntegral( q_initial );
   PrintDouble(res, "initial residual in NodeElimination");
   if(fabs(res) > tol) {
      bool SOL_FLAG = LeastSquaresNewton(ON, basis, q_temp, 0);
      if(SOL_FLAG == SOL_FOUND) {
         n_initial = q_temp->k;
         quadrature_realloc(q_temp->k, dim, dims, deg, q_new);
         quadrature_assign(q_temp, q_new);
      }
      else if(SOL_FLAG == SOL_NOT_FOUND) {
         Print("Initial quadrature did not converge. The initial guess should be more accurate.\n");
         FreeMemory(basis, q_temp, q_new);
         return;
      }
   }


   // run Node Elimination Algorithm. Theoretical optimum is reached when
   // (dim+1)*k = num_funcs, at which point elimination is pursued no further.
   int n_prev                   = -1;
   int n_cur                    = n_initial;
   bool_enum CONSTR_FLAG        = OFF;
   bool SOL_FLAG                = SOL_NOT_FOUND;
   double n_opt                 = ceil(1.0*num_funcs/(dim+1));
   double efficiency            = n_opt/n_initial;
   PrintElimInfo( dim, n_cur , n_opt, efficiency);
   while( ((dim+1)*n_cur > num_funcs)  && (n_cur >= 1) )
   {
      if(n_cur <= 1) goto FREERETURN;

#ifdef CONSTR_OPT
      int temp_prev;
      bool removed_node;
      if(dim == MAX_DIM)
         if(n_cur == n_prev) CONSTR_FLAG = ON;
         else                CONSTR_FLAG = OFF;
      else
         CONSTR_FLAG = OFF;
#endif
      n_prev = n_cur;

      // Initialize parameters and arrays
      quadrature_reinit(n_cur-1, q_temp);
      quadrature_reinit(n_cur, q_new);

      int nrows = num_funcs;
      int ncols = (dim+1) * n_cur;
      CMatrix J = CMatrix_init(nrows, ncols);
      GetJacobian(basis, q_new, J);
      CMatrix J_TR = CMatrix_init(J.cols, J.rows);
      CMatrix_Transpose(J, J_TR);
      // construct QR factorization of transpose of the JACOBIAN
      int INFO     = -1;
      int LDJ      = J_TR.rows;
      int N_REFL   = MIN(J_TR.rows, J_TR.cols);
      double *REFL = (double *)malloc( SIZE_DOUBLE(N_REFL) );
      double *WORK = (double *)malloc( SIZE_DOUBLE(J_TR.cols) );
      dgeqr2_(&J_TR.rows, &J_TR.cols, J_TR.id, &LDJ, REFL, WORK, &INFO);

      // initialize i*(dim+1)th rows of Q. Entries that correspond to weight indices are initialized to 1
      CMatrix QW  = CMatrix_init((dim+1)*n_cur, n_cur);
      for(int i = 0; i < QW.cols; ++i)
         QW.cid[i][i*(dim+1)] = 1.0;

      // obtain i*(dim+1)th rows of Q(from J_TR) explicitly and store them in QW columns
      DORMQR_M('L',  'T', REFL, J_TR, QW);
      COND_TEST_4;

      int *arrayIndex    = (int *)malloc( SIZE_INT(n_cur) );
      double *signif_ind = (double *)malloc( SIZE_DOUBLE(n_cur) );
      Vector Q2_ROW     = Vector_init(J_TR.rows);
      RMatrix Z          = RMatrix_init(n_cur, ncols);
      RMatrix dZ         = RMatrix_init(n_cur, ncols);
      // compute initial guesses for Newton's method and store them in Z
      for(i = 0; i < QW.cols; ++i)
      {
         // extract i*(dim+1)th row of Q2 and compute its norm, where Q = [Q1, Q2]
         for(j = 0; j < num_funcs; ++j) Q2_ROW.id[j] = 0.0;
         for(j = num_funcs; j < QW.rows; ++j) Q2_ROW.id[j] = QW.cid[i][j];
         double norm_Q2_ROW = TwoNorm(ncols-num_funcs, &Q2_ROW.id[num_funcs]);

         DORMQR_V('L', 'N', REFL, J_TR, Q2_ROW); // multiply Q by i*(dim+1)th row of Q2

         for(j = 0; j < ncols; ++j)
            dZ.rid[i][j] = Q2_ROW.id[j] * q_new->w[i]/SQUARE(norm_Q2_ROW);

         // compute new initial guesses and store them in predictor Z
         int count;
         for(count = 0, j = 0; j < n_cur; ++j)
         {
            Z.rid[i][(dim+1)*count] = q_new->w[j] - dZ.rid[i][(dim+1)*j];
            for(d = 0; d < dim; ++d)
               Z.rid[i][count*(dim+1)+1+d] = q_new->x[j*dim+d] - dZ.rid[i][j*(dim+1)+d+1];

            ++count;
         }

         signif_ind[i] = TwoNorm(ncols, dZ.rid[i]);
      }
      InsertionSort(n_cur, signif_ind, arrayIndex); // store indices of signif_ind in arrayIndex in ascending order

      free(signif_ind);
      free(REFL);
      free(WORK);
      Vector_free(Q2_ROW);
      RMatrix_free(dZ);
      CMatrix_free(J);
      CMatrix_free(J_TR);
      CMatrix_free(QW);


      // extract initial guesses and run Newton's method
      for(i = 0; i < n_cur; ++i)
      {
         int count;
         for(count = 0, j = 0; j < n_cur; ++j)
         {
            if(j == arrayIndex[i]) continue;
            for(d = 0; d < dim; ++d)
               q_temp->x[count*dim+d] = Z.rid[arrayIndex[i]][j*(dim+1) + d+1];

            q_temp->w[count] = Z.rid[arrayIndex[i]][j*(dim+1)];
            ++count;
         }

         if(V_InfNorm(q_temp->z) >= QUAD_HUGE) continue;

#ifdef CONSTR_OPT
         removed_node = false;
         ConstrVectData cVectData = ConstrVectDataInit();
         if(QuadInConstraint(q_temp) == false)
         {

            quadrature *q_prev = quadrature_init_full(n_cur-1, dim, dims, deg, D);
            for(count = 0, j = 0; j < n_cur; ++j)
            {
               if(j == arrayIndex[i]) continue;
               q_prev->w[count] = q_new->w[j];
               for(d = 0; d < dim; ++d)
                  q_prev->x[count*dim+d] = q_new->x[j*dim+d];
               ++count;
            }
            ConstrainedOptimization(q_prev, q_temp, &cVectData);

            if(cVectData.N_OR_W == WEIGHT)
            {
               temp_prev = q_temp->k;
               quadrature_remove_element(cVectData.boundaryNodeId, q_temp);
               removed_node = true;
            }

            quadrature_free(q_prev);
         }
#endif

         int its = 0;
         SOL_FLAG = LeastSquaresNewton(CONSTR_FLAG, basis, q_temp, &its);
         // store nodes and weights if Newton's method succeeded, update history
         if(SOL_FLAG == SOL_FOUND)
         {
            n_cur = q_temp->k;
            quadrature_realloc(q_temp->k, dim, dims, deg, q_new);
            quadrature_assign(q_temp, q_new);
            hist_data *hist_d = malloc(sizeof(hist_data));
            hist_d->nodes_tot = n_cur;
            hist_d->success_node = i;
            hist_d->success_its = its;
            glist_push(history, (void *)hist_d);
            break; // break i loop
         }
#ifdef CONSTR_OPT
         else if(SOL_FLAG == SOL_NOT_FOUND)
            if(removed_node == true)
               quadrature_realloc(temp_prev, dim, dims, deg, q_temp);
#endif
      }// end i loop
      RMatrix_free(Z);
      free(arrayIndex);

      if(SOL_FLAG == SOL_FOUND && CONSTR_FLAG == ON)
         Print("Solution found using constraint");

      // end NodeElimination if repeat flag is off and all iterations were unsuccessful
      if (SOL_FLAG == SOL_FOUND)
         PrintElimInfo(dim, n_cur, n_opt, n_opt/n_cur);
#ifdef CONSTR_OPT
      else if( (SOL_FLAG == SOL_NOT_FOUND) && (CONSTR_FLAG != ON)  && (dim == MAX_DIM) )
         Print("rerunning with constrained optimization");
#endif
      else if( SOL_FLAG == SOL_NOT_FOUND ) {
         Print("Last iteration did not converge");
         break;
      }

   }// end while loop

FREERETURN:
   // save nodes and weights
   q_final->k = n_cur;
   quadrature_assign(q_new, q_final);
   res = QuadTestIntegral(q_final);
   PrintDouble(res, "Final residual in NodeElimination"); printf("\n");

   FreeMemory(basis, q_temp, q_new);
}// end NodeElimination


static void FreeMemory(int_fast8_t *basis, quadrature *q_temp, quadrature *q_new)
{
   free(basis);
   quadrature_free(q_temp);
   quadrature_free(q_new);
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
   MAT_MUL(cols, rows, cols, Q_COPY_TR, Q_COPY, I);

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

