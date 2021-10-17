/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */

#include "NodeElimination.h"

#include "GetJacobian.h"
#include "InsertionSort.h"
#include "LeastSquaresNewton.h"
#include "TestIntegral.h"
#include "InDomain.h"
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

static void free_memory(int_fast8_t *basis, quadrature *q_temp, quadrature *q_new,
                        constraints *cons, ConstraintFuncs constrFuncs)
{
   free(basis);
   quadrature_free(q_temp);
   quadrature_free(q_new);
   constrFuncs.constr_free(cons);

}
static bool TestQR(int numRows, int qiCols, const double *Q);
static double TwoNorm(int n, double *z);

/* NodeElimination
 * Receives quadrature nodes, weights, domain and quadrature parameters and
 * uses the new node elimination scheme to eliminate one node to obtain initial guess for
 * Newton's method. Subsequently, routine calls Newton's method to obtain quadrature rule
 * with fewer nodes. The procedure is repeated until no more nodes can be eliminated.
 */
void NodeElimination(const quadrature *q_in, quadrature *q_final, const _DomainFuncs dom_funcs, const ConstraintFuncs constr_funcs, elim_history *hist)
{

   int n_initial = q_in->k;
   int num_funs = q_in->params->num_funs;
   int deg = q_in->params->deg;
   int dim = q_in->params->dim;
   int dim_plus_1 = dim+1;
   int *dims = q_in->params->dims;
   const DOMAIN_TYPE D = q_in->D;

   double tol = QUAD_TOL; // 10^(-15);
   double opt_factor = (double) num_funs/(dim_plus_1 * n_initial);

   quadrature *q_temp = quadrature_init(n_initial, dim, dims, deg, D);
   quadrature *q_new = quadrature_init(n_initial, dim, dims, deg, D);
   quadrature_assign(*q_in, *q_new);
   quadrature_assign(*q_in, *q_temp);

   int_fast8_t *basis = (int_fast8_t *)malloc( (num_funs*dim)*sizeof(int_fast8_t) );
   BasisIndices(deg, dim, basis);

   constraints *cons = constr_funcs.constr_init(dims);
   constr_funcs.get_constr(cons);

   // test accuracy of the initial quadrature
   // if residual is too large, attempt to find quadrature using Newton's method.
   // if Newton's method finds a solution, proceed to Node elimination, otherwise return.
   double res = TestIntegral( (const_quadrature *)q_in, dom_funcs );
   PRINT(res, "initial residual in NodeElimination");
   if( fabs(res) > tol )
   {
      int its = -1;
      BOOLEAN CONSTR_FLAG = OFF;
      bool SOL_FLAG = SOL_NOT_FOUND;
      SOL_FLAG = LeastSquaresNewton(CONSTR_FLAG, basis, dom_funcs, cons, &its, q_temp);
      if(SOL_FLAG == SOL_FOUND)
      {
         quadrature_assign(*q_temp, *q_new);
      }
      else if(SOL_FLAG == SOL_NOT_FOUND)
      {
         PRINT((char *)"Initial quadrature did not converge", 1);
         free_memory(basis, q_temp, q_new, cons, constr_funcs);
         return;
      }
   }

   // run Node Elimination Algorithm. Theoretical optimum is reached when
   // (dim+1)*k = num_funs, but a few more iterations are allowed in case under
   // exceptional circumstances theoretical optimum is surpassed.
   int n_prev = -1;
   int n_next = n_initial-1;
   int n_cur = n_initial;
   BOOLEAN REPEAT = OFF;
   BOOLEAN CONSTR_FLAG = OFF;
   bool SOL_FLAG = SOL_NOT_FOUND;
   while( (dim_plus_1*n_cur >= num_funs)  && (n_cur >= 1) )
   {
      if(n_cur <= 1) goto FREERETURN; // return if current number of nodes is 1

#ifdef CONSTR_OPT
      if(n_cur == n_prev)
      {
         CONSTR_FLAG = ON; REPEAT = ON;
      }
      else
      {
         CONSTR_FLAG = OFF; REPEAT = OFF;
      }
#endif
      n_prev = n_cur;

      // Initialize parameters and arrays
      quadrature_reinit(n_cur-1, q_temp);
      quadrature_reinit(n_cur, q_new);

      int ncols = dim_plus_1 * n_cur, nrows = num_funs, weight_cols = n_cur;

      Vector JACOBIAN = Vector_init(nrows*ncols);
      GetJacobian(dom_funcs, basis, *q_new, JACOBIAN.id);

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
      int nrows_Q = n_cur, ncols_Q = dim_plus_1 * n_cur;
      Vector QWEIGHT = Vector_init(nrows_Q*ncols_Q);
      memset(QWEIGHT.id, 0, SIZE_DOUBLE(nrows_Q*ncols_Q));
      for(int i = 0; i < nrows_Q; ++i)
         QWEIGHT.id[i*ncols + i*dim_plus_1] = 1.0;

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

            Z.id[i][dim_plus_1*count] = q_new->w[j] - dZ.id[i][dim_plus_1*j];
            for(int d = 0; d < dim; ++d)
               Z.id[i][count*dim_plus_1 + 1+d] = q_new->x[j*dim+d] - dZ.id[i][j*dim_plus_1 + d+1];

            ++count;
         }

         signif_ind[i] = TwoNorm(ncols, dZ.id[i]);
      }

      for(int j = 0; j < n_cur; ++j) arrayIndex[j] = j;
      InsertionSort(n_cur, signif_ind, arrayIndex); // store indices of signif_ind in arrayIndex in ascending order

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
               q_temp->x[j*dim+d] = Z.id[arrayIndex[i]][j*dim_plus_1 + d+1];
            q_temp->w[j] = Z.id[arrayIndex[i]][j*dim_plus_1];

         }

         int its = -1;
         SOL_FLAG = LeastSquaresNewton(CONSTR_FLAG, basis, dom_funcs, cons, &its, q_temp);
         // store nodes and weights if Newton's method succeeded, update history
         if(SOL_FLAG == SOL_FOUND)
         {
            n_cur = q_temp->k;
            q_new->k = n_cur;
            quadrature_assign(*q_temp, *q_new);
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
         PRINT((char *)"Solution found using constraint", 0);

      // end NodeElimination if repeat flag is off and all iterations were unsuccessful
      if (SOL_FLAG == SOL_FOUND)
         PrintElimInfo( dim, n_cur , (int)ceil(1.0*num_funs/(dim+1.0)) );
#ifdef CONSTR_OPT
      else if( (SOL_FLAG == SOL_NOT_FOUND) && (REPEAT != ON) )
      {
         PRINT((char *)"rerunning with constrained optimization", 3);
      }
#endif
      else if( SOL_FLAG == SOL_NOT_FOUND )
      {
         PRINT((char *)"Last iteration did not converge", 0);
         break;
      }

   }// end while loop

FREERETURN:
   // save nodes and weights
   q_final->k = n_cur;
   quadrature_assign(*q_new, *q_final);
   res = TestIntegral( (const_quadrature *)q_final, dom_funcs );
   PRINT(res, "Final residual in NodeElimination");
   printf("\n");

   free(basis);
   quadrature_free(q_temp);
   quadrature_free(q_new);
   constr_funcs.constr_free(cons);

}// end NodeElimination


static double TwoNorm(int n, double *z)
{
   assert(n >= 0);

   double norm = 0.0;
   for(int i = 0; i < n; ++i)
      norm += SQUARE(z[i]);

   norm = SQRT(norm);
   return norm;
}


static bool TestQR(int ncols, int qi_cols, const double *Q)
{
   char t1 = 'T', t2 = 'N';
   double alpha = 1.0, beta = 0.0;

   double max_non_diag = 0.0, max_diag = 1.0;
   double non_diag = 0.0, diag = 0.0;
   double temp = -1.0, tol = POW(10, -12);

   double *I = (double *)malloc( qi_cols*qi_cols*sizeof(double) );
   double *Q_COPY = (double *)malloc( ncols*qi_cols*sizeof(double) );
   memcpy(Q_COPY, Q, SIZE_DOUBLE(ncols*qi_cols));

   dgemm_(&t1, &t2, &qi_cols, &qi_cols, &ncols, &alpha, Q_COPY,
          &ncols, Q_COPY, &ncols, &beta, I, &qi_cols);

   for(int i = 0; i < qi_cols; ++i)
   {
      temp = fabs(I[i*qi_cols+i]);
      diag = MAX(temp, diag);

      for(int j = 0; j < qi_cols; ++j)
      {
         if(j == i) continue;
         temp = fabs(I[i*qi_cols+j]);
         non_diag = MAX(temp, non_diag);
      }

   }

   free(I);
   free(Q_COPY);

   double err1 = fabs(diag - max_diag);
   double err2 = fabs(non_diag - max_non_diag);
   double max_error = MAX(err1, err2);

   if(max_error > tol) return FAILED;
   else                return PASSED;
}

