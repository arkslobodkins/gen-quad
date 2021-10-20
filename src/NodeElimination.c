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


/* NodeElimination
 * Receives quadrature nodes, weights, domain and quadrature parameters and
 * uses the new node elimination scheme to eliminate one node to obtain initial guess for
 * Newton's method. Subsequently, routine calls Newton's method to obtain quadrature rule
 * with fewer nodes. The procedure is repeated until no more nodes can be eliminated.
 */
void NodeElimination(const quadrature *q_in, quadrature *q_final, elim_history *hist)
{
   if(q_in->setFuncsFlag != 1)
      PRINT_ERR("Error in NodeElimination, functions and constraints for "
                "q_in are not set", __LINE__, __FILE__);
   if(q_final->setFuncsFlag != 1)
      PRINT_ERR("Error in NodeElimination, functions and constraints for "
                "q_final are not set",  __LINE__, __FILE__);

   int n_initial       = q_in->k;
   int num_funs        = q_in->params->num_funs;
   int deg             = q_in->params->deg;
   int dim             = q_in->params->dim;
   int *dims           = q_in->params->dims;
   const DOMAIN_TYPE D = q_in->D;
   double tol          = QUAD_TOL; // 10^(-15);

   quadrature *q_temp = quadrature_init(n_initial, dim, dims, deg, D);
   quad_set_funcs_and_constr(q_temp);
   quadrature_assign(q_in, q_temp);

   quadrature *q_new = quadrature_init(n_initial, dim, dims, deg, D);
   quad_set_funcs_and_constr(q_new);
   quadrature_assign(q_in, q_new);

   int_fast8_t *basis = (int_fast8_t *)malloc( (num_funs*dim)*sizeof(int_fast8_t) );
   BasisIndices(deg, dim, basis);

   // test accuracy of the initial quadrature
   // if residual is too large, attempt to find quadrature using Newton's method.
   // if Newton's method finds a solution, proceed to Node elimination, otherwise return.
   double res = q_in->testIntegral( (const_quadrature *)q_in );
   PrintDouble(res, "initial residual in NodeElimination");
   if(fabs(res) > tol)
   {
      int its = -1; bool_enum CONSTR_FLAG = ON;
      bool SOL_FLAG = LeastSquaresNewton(CONSTR_FLAG, basis, q_temp, &its);
      if(SOL_FLAG == SOL_FOUND) {
         quadrature_assign(q_temp, q_new);
      }
      else if(SOL_FLAG == SOL_NOT_FOUND) {
         Print("Initial quadrature did not converge");
         FreeMemory(basis, q_temp, q_new);
         return;
      }
   }


   // run Node Elimination Algorithm. Theoretical optimum is reached when
   // (dim+1)*k = num_funs, but a few more iterations are allowed in case under
   // exceptional circumstances theoretical optimum is surpassed.
   int n_prev                   = -1;
   int n_cur                    = n_initial;
   bool_enum CONSTR_FLAG        = OFF;
   ATTR_UNUSED bool_enum REPEAT = OFF;
   bool SOL_FLAG                = SOL_NOT_FOUND;
   double n_opt                 = ceil(1.0*num_funs/(dim+1));
   double efficiency            = n_opt/n_initial;
   PrintElimInfo( dim, n_cur , n_opt, efficiency);
   while( ((dim+1)*n_cur >= num_funs)  && (n_cur >= 1) )
   {
      if(n_cur <= 1) goto FREERETURN; // return if current number of nodes is 1

#ifdef CONSTR_OPT
      if(dim == MAX_DIM) {
         if(n_cur == n_prev) {
            CONSTR_FLAG = ON; REPEAT = ON;
         }
         else {
            CONSTR_FLAG = OFF; REPEAT = OFF;
         }
      }
      else {
         CONSTR_FLAG = OFF; REPEAT = OFF;
      }

#endif
      n_prev = n_cur;

      // Initialize parameters and arrays
      quadrature_reinit(n_cur-1, q_temp);
      quadrature_reinit(n_cur, q_new);

      int ncols = (dim+1) * n_cur, nrows = num_funs, weight_cols = n_cur;

      Vector JACOBIAN = Vector_init(nrows*ncols);
      GetJacobian(basis, *q_new, JACOBIAN.id);

      // construct QR factorization of transpose of the JACOBIAN
      int INFO = -1;
      int J_trans_rows = ncols;
      int J_trans_cols = nrows;
      int LDJ = J_trans_rows;
      Vector *JACOBIAN_TRANS_F = &JACOBIAN; // Matrix in C is equivalent to its transpose in FORTRAN due to different memory layouts
      double *TAU = (double *) malloc( SIZE_DOUBLE(nrows) );
      double *WORK = (double *) malloc( SIZE_DOUBLE(ncols) );
      dgeqr2_(&J_trans_rows, &J_trans_cols, JACOBIAN_TRANS_F->id, &LDJ, TAU, WORK, &INFO);

      // initialize i*(dim+1)th rows of Q. Entries that correspond to weight indices are initialized to 1
      int nrows_Q = n_cur, ncols_Q = (dim+1) * n_cur;
      Vector QWEIGHT = Vector_init(nrows_Q*ncols_Q);
      memset(QWEIGHT.id, 0, SIZE_DOUBLE(nrows_Q*ncols_Q));
      for(int i = 0; i < nrows_Q; ++i)
         QWEIGHT.id[i*ncols + i*(dim+1)] = 1.0;

      // obtain i*(dim+1)th rows of Q explicitly and store them in QWEIGHT
      // QWEIGHT has c memory layout, since transpose was applied from the left side
      char TRANS = 'T'; char SIDE = 'L';
      int LDQ = J_trans_rows;
      int lworkQR = MAX(J_trans_rows, J_trans_cols);
      double *workQR = (double *) malloc( SIZE_DOUBLE(lworkQR));
      dormqr_(&SIDE, &TRANS, &J_trans_rows, &weight_cols, &J_trans_cols, JACOBIAN_TRANS_F->id, &LDJ, TAU, QWEIGHT.id, &LDQ, workQR, &lworkQR, &INFO);
      COND_TEST_4;



      Vector Q2_W_ROW = Vector_init(ncols);
      int *arrayIndex = (int *) malloc( SIZE_INT(n_cur) );
      double *signif_ind = (double *) malloc( SIZE_DOUBLE(n_cur) );
      Matrix Z = Matrix_init(n_cur, ncols);
      Matrix dZ = Matrix_init(n_cur, ncols);
      // compute initial guesses for Newton's method and store them in Z
      for(int i = 0; i < nrows_Q; ++i)
      {
         // extract i*(dim+1)th row of Q2 and compute its norm, where Q = [Q1, Q2]
         for(int j = 0; j < num_funs; ++j) Q2_W_ROW.id[j] = 0.0;
         for(int j = num_funs; j < ncols_Q; ++j) Q2_W_ROW.id[j] = QWEIGHT.id[i*ncols+j];
         double norm_Q2_W_ROW = TwoNorm(ncols-num_funs, &Q2_W_ROW.id[num_funs]);

         // multiply Q by i*(dim+1)th row of Q2
         TRANS = 'N'; char SIDE = 'L';
         int ONE_COLUMN = 1;
         dormqr_(&SIDE, &TRANS, &J_trans_rows, &ONE_COLUMN, &J_trans_cols, JACOBIAN_TRANS_F->id, &LDJ, TAU, Q2_W_ROW.id, &LDQ, workQR, &lworkQR, &INFO);

         for(int j = 0; j < ncols; ++j)
         {
            dZ.id[i][j] = Q2_W_ROW.id[j];
            dZ.id[i][j] = dZ.id[i][j] * q_new->w[i]/SQUARE(norm_Q2_W_ROW);
         }

         // compute new initial guesses and store them in predictor Z
         int count = 0;
         for(int j = 0; j < n_cur; ++j)
         {
            if(i == j) continue;

            Z.id[i][(dim+1)*count] = q_new->w[j] - dZ.id[i][(dim+1)*j];
            for(int d = 0; d < dim; ++d)
               Z.id[i][count*(dim+1) + 1+d] = q_new->x[j*dim+d] - dZ.id[i][j*(dim+1) + d+1];

            ++count;
         }

         signif_ind[i] = TwoNorm(ncols, dZ.id[i]);
      }
      InsertionSort(n_cur, signif_ind, arrayIndex); // store indices of signif_ind in arrayIndex in ascending order
//      for(int s = 0; s < n_cur; ++s)
//         arrayIndex[s] = s;

      free(signif_ind);
      free(TAU);
      free(WORK);
      free(workQR);
      Matrix_free(dZ);
      Vector_free(Q2_W_ROW);
      Vector_free(JACOBIAN);
      Vector_free(QWEIGHT);

      // extract initial guesses and run Newton's method
      for(int i = 0; i < n_cur; ++i)
      {
         for(int j = 0; j < n_cur-1; ++j)
         {
            for(int d = 0; d < dim; ++d)
               q_temp->x[j*dim+d] = Z.id[arrayIndex[i]][j*(dim+1) + d+1];

            q_temp->w[j] = Z.id[arrayIndex[i]][j*(dim+1)];
         }

         int its = -1;

         ConstrNodeData cNodeData;
         if(q_temp->inConstraint((const_quadrature *)q_temp) == false)
         {
            quadrature *q_prev = quadrature_init(n_cur-1, dim, dims, deg, D);
            int count = 0;
            for(int r = 0; r < n_cur; ++r)
            {
               if(r == arrayIndex[i]) continue;
               q_prev->w[count] = q_new->w[r];
               for(int d = 0; d < dim; ++d)
                  q_prev->x[count*dim+d] = q_new->x[r*dim+d];
               ++count;
            }

            cNodeData = constrain_vector(q_temp->FULL_A, q_temp->FULL_b,
                  (const_quadrature *)q_prev, (const_quadrature *)q_temp);
            Vector q_diff = Vector_init((n_cur-1)*(dim+1));
            for(int l = 0; l < q_diff.len; ++l)
               q_diff.id[l] = q_prev->z[l] - q_temp->z[l];
            for(int l = 0; l < q_diff.len; ++l)
               q_temp->z[l] = q_prev->z[l] + (-cNodeData.tMin + POW(10, -12))*q_diff.id[l];


            Vector_free(q_diff);
            quadrature_free(q_prev);
         }


         SOL_FLAG = LeastSquaresNewton(CONSTR_FLAG, basis, q_temp, &its);
         // store nodes and weights if Newton's method succeeded, update history
         if(SOL_FLAG == SOL_FOUND)
         {
            n_cur = q_temp->k;
            q_new->k = n_cur;
            quadrature_assign(q_temp, q_new);
            hist->nodes_tot[hist->tot_elims] = n_cur;
            hist->success_node[hist->tot_elims] = i;
            hist->success_its[hist->tot_elims] = its;
            ++hist->tot_elims;
            break; // break i loop
         }

      }// end i loop
      Matrix_free(Z);
      free(arrayIndex);



      if(SOL_FLAG == SOL_FOUND && CONSTR_FLAG == ON)
         Print("Solution found using constraint");

      // end NodeElimination if repeat flag is off and all iterations were unsuccessful
      if (SOL_FLAG == SOL_FOUND)
         PrintElimInfo(dim, n_cur, n_opt, n_opt/n_cur);
#ifdef CONSTR_OPT
      else if( (SOL_FLAG == SOL_NOT_FOUND) && (REPEAT != ON)  && (dim == MAX_DIM) )
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
   res = q_final->testIntegral( (const_quadrature *)q_final );
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
   assert(n >= 0);

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
   double temp = -1.0, tol = POW(10, -12);

   double *Q_COPY = (double *)calloc( rows*cols, sizeof(double) );
   memcpy(Q_COPY, Q, SIZE_DOUBLE(rows*cols));

   double *Q_COPY_TRANS = (double *)calloc( cols*rows, sizeof(double) );
   Transpose(cols, rows, Q_COPY, Q_COPY_TRANS);

   double *I = (double *)calloc( cols*cols, sizeof(double) );
   MAT_MUL(cols, rows, cols, Q_COPY_TRANS, Q_COPY, I);

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
   free(Q_COPY_TRANS);

   double err1 = fabs(diag - max_diag);
   double err2 = fabs(non_diag - max_non_diag);
   double max_error = MAX(err1, err2);

   if(max_error > tol) return FAILED;
   else                return PASSED;
}

