/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */

// for timing LAPACK routine globally, some profiling tools lack this info
#include <time.h>
double LSQ_TIME = 0.0;

#include "LeastSquaresNewton.h"

#include "GetFunction.h"
#include "GetJacobian.h"
#include "TestIntegral.h"
#include "Constraints.h"
#include "Vector.h"
#include "Matrix.h"
#include "Quadrature.h"
#include "GENERAL_QUADRATURE.h"
#include "LINALG.h"
#include "Print.h"
#include "Conditional_Debug.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>

#define MAX_ELIM_WEIGHTS 10


//static inline bool InConstraintElem(const _DomainFuncs dom_funcs, const_quadrature *quad, int elem);
//static inline bool PosWeightElem(const_quadrature *q, int elem);
static bool PosWeights(const_quadrature *q);
static bool InConstraint(const _DomainFuncs funcs, const_quadrature *quad);
static bool CheckForFail(int info, double error_norm, double error_norm_prev, Vector least_sq_sol);
static void Transpose(int M, int N, const double *A, double *B);


#ifdef CONSTR_OPT
static constr_vect_data constr_vect_data_init(int n_cols);
static void constr_vect_data_reset(int n_cols, constr_vect_data  *C_V_DATA);
// ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  || || ||  ||  ||  ||  ||  ||  ||  ||  || ||
//static void constr_vect_data_realloc(int n_cols, constr_vect_data *C_V_DATA);
//static void constr_vect_data_free(constr_vect_data C_V_DATA);
//static void AppendJacobian(const int nRows, const int nCols, const double *eqnConstr, Vector jacobian);
//static void AppendFunction(const int nRows, double *f);
// ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  || || ||  ||  ||  ||  ||  ||  ||  ||  || ||


static void ConstrainedProjection(const _DomainFuncs dom_funcs, const_quadrature *q_prev, quadrature *q_next, const constraints *cons);
static bool ConstrainedOptimization(const _DomainFuncs dom_funcs, const constraints *cons, const_quadrature *q_prev, quadrature *q_next, constr_vect_data *C_V_DATA);

static constr_node_data constrain_vector(const Matrix A, const Vector b, const_quadrature *q_prev, const_quadrature *q_next);


static constr_node_data ShortenNode(const Matrix A, const Vector b, const Vector z_old, const Vector z_new);
static void ProjectNode(const Matrix eqn_matrix, const Vector dx, Vector x_projected);
#endif




/* LeastSquaresNewton
 * Receives initial quadrature guess. Primarily solves
 * underdetermined systems of equations in the least
 * squares sense. Returns success if algorithm converged
 * and all nodes are inside of the domain and if all nodes are positive..
 */
bool LeastSquaresNewton(const BOOLEAN FLAG_CONSTR, const int_fast8_t *basis, const _DomainFuncs dom_funcs,
      const constraints *cons, int *its, quadrature *q_orig)
{
   assert(q_orig->k >= 0);
   int k = q_orig->k;

   int WEIGHT_ELIM_FLAG = 0;
   bool SOL_FLAG = SOL_NOT_FOUND;
   int p = q_orig->params->deg;
   int num_funcs = q_orig->params->num_funs;
   int dim = q_orig->params->dim;
   int *dims = q_orig->params->dims;
   const DOMAIN_TYPE D = q_orig->D;

   // quadrature nodes and weights expressed as vectors, more convenient
   // to use for vector operations than quadrature objects
   int n_rows = num_funcs, n_cols = (dim+1)*k;
   Vector jacobian = Vector_init(n_rows * n_cols);
   Vector jacobian_transpose = Vector_init(n_rows * n_cols);
   Vector rhs = Vector_init(num_funcs);

   quadrature *q_prev = quadrature_init(k, dim, dims, p, D);
   quadrature *q_next = quadrature_init(k, dim, dims, p, D);
   quadrature_assign(*q_orig, *q_prev);
   quadrature_assign(*q_orig, *q_next);

   int its_loc = 0; int maxiter = 25;
   double error_norm = -1.0, error_norm_prev = -1.0, error_norm_update = -1;
   double tol = QUAD_TOL; // 10^(-15);


   // initialize LAPACK Parameters
   int info = 0, nRhs = 1, lda = n_rows, lead_dim = MAX(n_rows, n_cols),
       lwork = 5*n_cols, small_dim = MIN(n_rows, n_cols);
   char trans = 'N';
   Vector work = Vector_init(lwork);
   Vector least_sq_sol = Vector_init(lead_dim);


#ifdef CONSTR_OPT

   int elim_weights = 0;
   bool CONSTR_RETURN = CONSTRAINT_FAILURE;
   constr_vect_data C_V_DATA = constr_vect_data_init(n_cols);
#else
   if(FLAG_CONSTR == ON)
      PRINT_ERR("constrained optimization is turned off and won't be performed", 1, __LINE__, __FILE__);

#endif


   // return if input is a satisfactory quadrature
   GetFunction(dom_funcs, basis, *q_prev, rhs.id);
   error_norm = V_ScaledTwoNorm(rhs);
   if( (error_norm < tol) && (InConstraint( dom_funcs, (const_quadrature *)q_prev ) == true) )
   {
      SOL_FLAG = SOL_FOUND;
      *its = 0;
      goto FREERETURN;
   }


   do
   {

#ifdef CONSTR_OPT
      switch(C_V_DATA.ACTIVE_CONSTRAINTS)
      {

         case OFF:
            n_rows = num_funcs;
            n_cols = (dim+1)*k;
            break;

         case ON:
            switch(C_V_DATA.N_OR_W)
            {

               case NODE:
                  n_rows = num_funcs;
                  n_cols = (dim+1)*k;
// ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  || || ||  ||  ||  ||  ||  ||  ||  ||  || ||
//                  n_rows = num_funcs + 1;
//                  n_cols = (dim+1)*k;
// ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  || || ||  ||  ||  ||  ||  ||  ||  ||  || ||
                  break;

               case WEIGHT:

                  ++elim_weights;
                  if(k == 1  || elim_weights > MAX_ELIM_WEIGHTS)
                  {
                     SOL_FLAG = SOL_NOT_FOUND;
                     goto FREERETURN;
                  }

                  WEIGHT_ELIM_FLAG = 1;
                  --k;
                  n_cols = (dim+1)*k;
                  n_rows = num_funcs;
                  quadrature_remove_element(C_V_DATA.boundary_node_id, q_next);
                  quadrature_remove_element(C_V_DATA.boundary_node_id, q_prev);
// ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  || || ||  ||  ||  ||  ||  ||  ||  ||  || ||
//                  constr_vect_data_realloc((dim+1)*k, &C_V_DATA);
// ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  || || ||  ||  ||  ||  ||  ||  ||  ||  || ||
                  break;

               default:
                  break;

            }
      }
#endif

      // update LAPACK parameters. Their size varies depending on whether ACTIVE_CONSTRAINTS flag was set or not
      // during the previous iteration.
      small_dim = MIN(n_rows, n_cols);
      lda = n_rows;
      lead_dim = MAX(n_rows, n_cols);
      Vector_realloc(n_rows*n_cols, &jacobian);
      Vector_realloc(n_rows*n_cols, &jacobian_transpose);
      Vector_realloc(n_rows, &rhs);
      Vector_realloc(lead_dim, &least_sq_sol);

      // Compute function and jacobian of a function to be solved
      GetJacobian(dom_funcs, basis, *q_prev, jacobian.id);
      GetFunction(dom_funcs, basis, *q_prev, rhs.id);

// ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  || || ||  ||  ||  ||  ||  ||  ||  ||  || ||
//#ifdef CONSTR_OPT
//      if ( (C_V_DATA.ACTIVE_CONSTRAINTS == ON )  && (C_V_DATA.N_OR_W == NODE) ) // add additional equation if constrained optimization is set
//      {
//         AppendJacobian(n_rows, n_cols, C_V_DATA.additional_constr_eqn.id, jacobian);
//         AppendFunction(n_rows, rhs.id);
//      }
//#endif
// ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  || || ||  ||  ||  ||  ||  ||  ||  ||  || ||

      for(int i = 0; i < small_dim; ++i)
         least_sq_sol.id[i] = rhs.id[i];
      for(int i = small_dim; i < lead_dim; ++i)
         least_sq_sol.id[i] = 0.0;

      // convert memory layout of jacobian to that of LAPACK
      C_TO_FORTRAN(n_rows, n_cols, jacobian.id, jacobian_transpose.id);
      // solve least squares problem using LAPACK routine and time time it
      time_t start = clock();
      dgels_(&trans, &n_rows, &n_cols, &nRhs, &jacobian_transpose.id[0], &lda,
            &least_sq_sol.id[0], &lead_dim, &work.id[0], &lwork, &info);
      time_t end = clock();
      LSQ_TIME += (double)(end - start)/CLOCKS_PER_SEC;

      for(int i = 0; i < n_cols; ++i)
         q_next->z[i] = q_prev->z[i] - least_sq_sol.id[i];

      error_norm_prev = error_norm;
      GetFunction(dom_funcs, basis, *q_next, rhs.id);
      error_norm = V_ScaledTwoNorm(rhs);

      bool check_values = CheckForFail(info, error_norm, error_norm_prev, least_sq_sol);
      if(check_values == true)
      {
         SOL_FLAG = SOL_NOT_FOUND;
         goto FREERETURN;
      }



#ifdef CONSTR_OPT
      if(FLAG_CONSTR == ON)
      {
         ConstrainedProjection(dom_funcs, (const_quadrature *)q_prev, q_next, cons);
         CONSTR_RETURN = ConstrainedOptimization(dom_funcs, cons, (const_quadrature *)q_prev, q_next, &C_V_DATA);

         if(CONSTR_RETURN == CONSTRAINT_FAILURE)
         {
            SOL_FLAG = SOL_NOT_FOUND;
            goto FREERETURN;
         }
      }
#endif


      GetFunction(dom_funcs, basis, *q_next, rhs.id);
      error_norm_update = V_ScaledTwoNorm(rhs);
      quadrature_assign(*q_next, *q_prev);
      ++its_loc;

   } while( (its_loc < maxiter) && (error_norm_update > tol) );


   // check if quadrature satisfies constraints
   if( (InConstraint( dom_funcs,  (const_quadrature *)q_next ) == false) || (error_norm_update > tol) )
   {
      SOL_FLAG = SOL_NOT_FOUND;
      goto FREERETURN;
   }

   // check whether Newton's method succeeded
   if( (its_loc < maxiter)
         && (isnan(error_norm_update) == 0)
         && (isinf(error_norm_update) == 0)
         && (error_norm_update <= tol) )
   {
      q_orig->k = k;
      quadrature_assign(*q_next, *q_orig);
      *its = its_loc;
      SOL_FLAG = SOL_FOUND;
      if(WEIGHT_ELIM_FLAG == 1)
         PRINT((char *)"SUCCEEDED AFTER ELIMINATING WEIGHT", 0);
   }
   else
   {
      *its = maxiter;
      SOL_FLAG = SOL_NOT_FOUND;
   }

FREERETURN:
   Vector_free(least_sq_sol);
   Vector_free(work);
   quadrature_free(q_prev);
   quadrature_free(q_next);
   Vector_free(rhs);
   Vector_free(jacobian_transpose);
   Vector_free(jacobian);

// ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  || || ||  ||  ||  ||  ||  ||  ||  ||  || ||
#ifdef CONSTR_OPT
//   constr_vect_data_free(C_V_DATA);
#endif
// ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  || || ||  ||  ||  ||  ||  ||  ||  ||  || ||

   return SOL_FLAG;
}// end LeastSquaresNewton


static bool InConstraint(const _DomainFuncs dom_funcs, const_quadrature *quad)
{
   return dom_funcs.inDomain(quad) & PosWeights(quad);
}


//static inline bool InConstraintElem(const domainFuncs dom_funcs, const_quadrature *quad, int elem)
//{
//   return dom_funcs.inDomainElem(quad, elem) & PosWeightElem(quad, elem);
//}
//
//
//static inline bool PosWeightElem(const_quadrature *q, int elem)
//{
//   return q->w[elem] > 0 ? true : false;
//}


static bool PosWeights(const_quadrature *q)
{

   for(int i = 0; i < q->k; ++i)
      if(q->w[i] < 0)
         return false;

   return true;
}


static bool CheckForFail(int info, double error_norm, double error_norm_prev, Vector least_sq_sol)
{

   if( error_norm > (error_norm_prev+1) )      // return if Newton's method does not have reasonable convergence
      return 1;

   else if( V_InfNorm(least_sq_sol) > 20 )    // return if solution is exploding
      return 1;

   else if(info > 0)                           // return if LAPACK routine has failed
      return 1;

   else if( V_CheckInfAndNan(least_sq_sol) )
      return 1;

   else return 0;

}


static void Transpose(int M, int N, const double *A, double *B)
{
   for(int i = 0; i < M; ++i)
   {
      for(int j = 0; j < N; ++j)
         B[i+j*M] = A[j+i*N];
   }
}




#ifdef CONSTR_OPT

static constr_vect_data constr_vect_data_init(int n_cols)
{
   constr_vect_data C_V_DATA;
   C_V_DATA.ACTIVE_CONSTRAINTS = OFF;
   C_V_DATA.N_OR_W = NONE;
   C_V_DATA.boundary_node_id = -1;
   C_V_DATA.eqn_id = -1;


// ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  || || ||  ||  ||  ||  ||  ||  ||  ||  || ||
//   Vector_init(n_cols, &C_V_DATA->additional_constr_eqn);
// ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  || || ||  ||  ||  ||  ||  ||  ||  ||  || ||
   return C_V_DATA;
}


static void constr_vect_data_reset(int n_cols, constr_vect_data *C_V_DATA)
{
   C_V_DATA->ACTIVE_CONSTRAINTS = OFF;
   C_V_DATA->N_OR_W = NONE;
   C_V_DATA->boundary_node_id = -1;
   C_V_DATA->eqn_id = -1;

// ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  || || ||  ||  ||  ||  ||  ||  ||  ||  || ||
//   memset(C_V_DATA->additional_constr_eqn.id, 0, n_cols*sizeof(double));
// ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  || || ||  ||  ||  ||  ||  ||  ||  ||  || ||
}


// ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  || || ||  ||  ||  ||  ||  ||  ||  ||  || ||
//static void constr_vect_data_realloc(int n_cols, constr_vect_data *C_V_DATA)
//{
//
//   Vector_realloc(n_cols, &C_V_DATA->additional_constr_eqn);
//}
// ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  || || ||  ||  ||  ||  ||  ||  ||  ||  || ||


// ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  || || ||  ||  ||  ||  ||  ||  ||  ||  || ||
//static void constr_vect_data_free(constr_vect_data C_V_DATA)
//{
//   Vector_free(C_V_DATA.additional_constr_eqn);
//}
// ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  || || ||  ||  ||  ||  ||  ||  ||  ||  || ||


// ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  || || ||  ||  ||  ||  ||  ||  ||  ||  || ||
//static void AppendJacobian(int n_rows, int n_cols, const double *eqnConstr, Vector jacobian)
//{
//   assert(jacobian.len == n_rows*n_cols);
//   int last_row = n_rows-1;
//   for(int j = 0; j < n_cols; ++j)
//      jacobian.id[ij2(last_row, j, n_cols)] = eqnConstr[j];
//
//}
//
//
//static void AppendFunction(int len, double *f)
//{
//   f[len-1] = 0.0;
//}
// ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  || || ||  ||  ||  ||  ||  ||  ||  ||  || ||


static void ConstrainedProjection(const _DomainFuncs dom_funcs, const_quadrature *q_prev, quadrature *q_next, const constraints *cons)

{
   int k = q_prev->k;
   int dim = q_prev->params->dim;
   assert(dim == cons->M.cols);
   double tol = pow(10, -12);

   int num_eqns = cons->M.rows;
   Matrix M = cons->M;
   Vector b = cons->b;

   Vector node_change = Vector_init(dim);
   Vector node_projected = Vector_init(dim);

   for(int i = 0; i < k; ++i)
   {
      int active_eqn_count = 0, node_index = i*dim;
      bool do_project = false;
      bool eqn_flags[num_eqns]; for(int n = 0; n < num_eqns; ++n) eqn_flags[n] = false;

      for(int j = 0; j < num_eqns; ++j)
      {
         double lhs = 0.0, lhs_prev = 0.0;
         for(int d = 0; d < dim; ++d)
            lhs_prev += M.id[j][d] * q_prev->x[node_index+d];
         for(int d = 0; d < dim; ++d)
            lhs += M.id[j][d] * q_next->x[node_index+d];

         if( fabs( lhs_prev - b.id[j] ) <= tol )
         {
            ++active_eqn_count;
            eqn_flags[j] = true;
            if( ! do_project)
               if( lhs  >= b.id[j] ) do_project = true;
         }
      }

      if(active_eqn_count > 0  && do_project == true)
      {
         int count = 0;
         Matrix eqn_matrix = Matrix_init(dim, active_eqn_count);

         for(int j = 0; j < num_eqns; ++j)
         {
            if(eqn_flags[j] == true)
            {
               for(int d = 0; d < dim; ++d)
                  eqn_matrix.id[d][count]= M.id[j][d];

               ++count;
            }
         }

         for(int s = 0; s < dim; ++s)
            node_change.id[s] = q_next->x[node_index+s] - q_prev->x[node_index+s];

         ProjectNode(eqn_matrix, node_change, node_projected);
         for(int s = 0; s < dim; ++s)
            q_next->x[node_index+s] = q_prev->x[node_index+s] + node_projected.id[s];

         Matrix_free(eqn_matrix);
      }

   }
   COND_TEST_1;

   Vector_free(node_change);
   Vector_free(node_projected);
}


static bool ConstrainedOptimization(const _DomainFuncs dom_funcs, const constraints *cons, const_quadrature *q_prev, quadrature *q_next, constr_vect_data *C_V_DATA)
{
   assert(q_prev->k == q_next->k);

   int k = q_prev->k;
   int dim = q_prev->params->dim;
   int n_cols = (dim+1)*k;


   bool constrNext_id[k], constr_prev_id[k];

   double q_diff_id[n_cols];
   Vector q_diff = {0};
   q_diff.len = k; q_diff.id = q_diff_id;
   for(int i = 0; i < n_cols; ++i)
      q_diff.id[i] = q_prev->z[i] - q_next->z[i];

   int A_rows = cons->M.rows+1, A_cols = cons->M.cols+1;
   Matrix A = {0};
   A.rows = A_rows; A.cols = A_cols;
   double A_id[A_rows][A_cols]; double *A_id_ptr[A_rows]; A.id = A_id_ptr;
   for(int i = 0; i < A.rows; ++i)
   {
      A.id[i] = A_id[i];
      memset(A.id[i], 0, A_cols*sizeof(double));
   }

   for(int i = 0; i < A.rows-1; ++i)
      for(int j = 0; j < A.cols-1; ++j)
         A.id[i+1][j+1] = cons->M.id[i][j];
   A.id[0][0] = -1.0;

   Vector b = {0};
   b.len = cons->b.len+1;
   double b_id[b.len]; b.id = b_id;
   for(int i = 0; i < b.len-1; ++i)
      b.id[i+1] = cons->b.id[i];
   b.id[0] = 0.0;

   // Compute only if equations at the previous iteration satisfy constraints
   if( ( InConstraint(dom_funcs, (const_quadrature *)q_next) == false )
         && ( InConstraint(dom_funcs, q_prev ) == true ) )
   {
      constr_node_data constr_node_d = constrain_vector(A, b, q_prev, (const_quadrature *)q_next);
      COND_TEST_2;

      if(constr_node_d.flag == 0)
         return CONSTRAINT_FAILURE;

      switch(constr_node_d.N_OR_W)
      {

         case NODE:
            for(int i = 0; i < n_cols; ++i)
               q_next->z[i] = q_prev->z[i] + (-constr_node_d.t_min + POW(10, -12) ) * q_diff.id[i];

// ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  || || ||  ||  ||  ||  ||  ||  ||  ||  || ||
//            // compute entries of additional equation for constraining the boundary
//            int index = -1;
//            int eqn = constr_node_d.eqn_id;
//            memset(C_V_DATA->additional_constr_eqn.id, 0, n_cols * sizeof(double));
//            for(int j = 0; j < cons.M.cols; ++j)
//            {
//               if(j == 0)  index = constr_node_d.node_id;
//               else        index = k + constr_node_d.node_id * dim-1;
//               C_V_DATA->additional_constr_eqn.id[index+j] = cons.M.id[eqn][j];
//            }
// ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  || || ||  ||  ||  ||  ||  ||  ||  ||  || ||

            C_V_DATA->ACTIVE_CONSTRAINTS = ON;
            C_V_DATA->N_OR_W = NODE;
            C_V_DATA->boundary_node_id = constr_node_d.node_id;
            C_V_DATA->eqn_id = constr_node_d.eqn_id;

            break;

         case WEIGHT:
            for(int i = 0; i < n_cols; ++i)
               q_next->z[i] = q_prev->z[i] + (-constr_node_d.t_min + POW(10, -13) ) * q_diff.id[i];

// ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  || || ||  ||  ||  ||  ||  ||  ||  ||  || ||
//            memset(C_V_DATA->additional_constr_eqn.id, 0, n_cols*sizeof(double));
// ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  ||  || || ||  ||  ||  ||  ||  ||  ||  ||  || ||
            C_V_DATA->ACTIVE_CONSTRAINTS = ON;
            C_V_DATA->N_OR_W = WEIGHT;
            C_V_DATA->boundary_node_id = constr_node_d.node_id;
            C_V_DATA->eqn_id = constr_node_d.eqn_id;
            break;

         default:
            constr_vect_data_reset(n_cols, C_V_DATA);
            break;
      }
      COND_TEST_3;

      assert(InConstraint(dom_funcs, (const_quadrature *)q_next) == 1);

   }
   else
      constr_vect_data_reset(n_cols, C_V_DATA);



   return CONSTRAINT_SUCCESS;
}


static constr_node_data constrain_vector(const Matrix A, const Vector b, const_quadrature *q_prev, const_quadrature *q_next)
{
   assert(A.cols == q_prev->params->dim+1);
   assert(q_prev->k == q_next->k);

   int k = q_prev->k;
   int dim = q_prev->params->dim;
   constr_node_data *C_N_DATA_loc = (constr_node_data *)malloc( k* sizeof(constr_node_data) );

   int *eqn = (int *)calloc( k, sizeof(int) );
   for(int i = 0; i < k; ++i)
   {
      C_N_DATA_loc[i].N_OR_W = NONE;
      C_N_DATA_loc[i].eqn_id = -1;
      C_N_DATA_loc[i].node_id = -1;
      C_N_DATA_loc[i].flag = 0;
      C_N_DATA_loc[i].t_min = 1.0;
   }

   Vector z_prev_node = {0}, z_next_node = {0};
   z_prev_node.len = dim+1; z_next_node.len = dim+1;
   double z_prev_id[dim+1], z_next_id[dim+1];
   z_prev_node.id = z_prev_id; z_next_node.id = z_next_id;

   int i, count;
   for(count = 0, i = 0; i < k; ++i)
   {
      z_prev_node.id[0] = q_prev->w[i];
      z_next_node.id[0] = q_next->w[i];
      for(int d = 0; d < dim; ++d)
      {
         z_prev_node.id[d+1] = q_prev->x[i*dim+d];
         z_next_node.id[d+1] = q_next->x[i*dim+d];
      }

      double lhs[A.rows];
      for(int r = 0; r < A.rows; ++r)
         lhs[r] = A.id[r][0] * q_next->w[i];

      for(int r = 0; r < A.rows; ++r)
         for(int c = 0; c < A.cols-1; ++c)
            lhs[r] += A.id[r][c+1] * q_next->x[i*dim+c];

      bool shorten = false;
      for(int r = 0; r < A.rows; ++r)
      {
         if(lhs[r] > b.id[r])
         {
            shorten = true;
            break;
         }
      }


      if( shorten )
      {
         C_N_DATA_loc[i] = ShortenNode(A, b, z_prev_node, z_next_node);

         if( C_N_DATA_loc[i].flag == 1 ) // store equations if succeeded
            ++count;
      }
   }


   Vector t_vec = Vector_init(k);
   for(int i = 0; i < k; ++i)
      t_vec.id[i] = C_N_DATA_loc[i].t_min;

   int min_index = -1;
   V_min v_min = {0};
   constr_node_data C_N_DATA;

   if(count > 0) // compute minimum t for mapping nodes to the boundary, defaults to 1
   {
      v_min = VectorMin(t_vec);
      min_index = v_min.min_index;
      C_N_DATA.N_OR_W = C_N_DATA_loc[min_index].N_OR_W;
      C_N_DATA.eqn_id = C_N_DATA_loc[min_index].eqn_id;
      C_N_DATA.node_id = min_index;
      C_N_DATA.t_min = v_min.min_value;
      C_N_DATA.flag = 1;

      assert(v_min.min_value < 1.0001  && v_min.min_value >= 0);
   }
   else
   {
      C_N_DATA.N_OR_W = NONE;
      C_N_DATA.eqn_id = -1;
      C_N_DATA.node_id = -1;
      C_N_DATA.t_min = 1.0;
      C_N_DATA.flag = 0;

   }

   free(C_N_DATA_loc);
   free(eqn);
   Vector_free(t_vec);

   return C_N_DATA;
}


static constr_node_data ShortenNode(const Matrix A, const Vector b_bound, const Vector z_old, const Vector z_new)
{
   assert(A.rows == b_bound.len);

   int i = -1, j = -1;
   int out_count = -1;
   int n_rows = A.rows;
   int n_cols = A.cols;
   int coord_num[n_rows];
   constr_node_data cnd;

   // Allocate vectors on the stack
   Vector b_old, b_new, b_diff, t, dz;
   b_old.len = n_rows; b_new.len = n_rows; b_diff.len = n_rows; t.len = n_rows; dz.len = n_cols;
   double b_old_id[n_rows], b_new_id[n_rows], b_diff_id[n_rows], t_id[n_rows], dz_id[n_cols];
   b_old.id = b_old_id; b_new.id = b_new_id; b_diff.id = b_diff_id; t.id = t_id; dz.id = dz_id;
   memset(t.id, 0, n_rows*sizeof(double));
   memset(b_old.id, 0, n_rows*sizeof(double));
   memset(b_diff.id, 0, n_rows*sizeof(double));
   memset(b_new.id, 0, n_rows*sizeof(double));
   memset(coord_num, -1, n_rows*sizeof(int));


   V_AddScale(1.0, z_new, -1.0, z_old, dz);
   MatVec(A, z_new, b_new);

   for (out_count = 0, i = 0; i < n_rows; ++i)
   {
      if (b_new.id[i] > b_bound.id[i])
      {
         for (j = 0; j < n_cols; ++j)
            b_old.id[i] += A.id[i][j] * z_old.id[j];
         for (j = 0; j < n_cols; ++j)
            b_diff.id[i] += A.id[i][j] * dz.id[j];

         t.id[i] = (b_bound.id[i] - b_old.id[i]) / b_diff.id[i];
         if(t.id[i] > 0.0)
         {
            coord_num[out_count] = i;
            ++out_count;
         }
      }
   }

   // exit if all nodes already satisfy constraints
   if(out_count == 0)
   {
      cnd.eqn_id = -1;
      cnd.node_id = -1;
      cnd.flag = 0;
      cnd.t_min = 1.0;
      cnd.N_OR_W = NONE;
      return cnd;
   }
   else
   {
      // compute minimum t
      int eqn_id = coord_num[0];
      double t_min = t.id[coord_num[0]];
      for(i = 1; i < out_count; ++i)
      {
         if (t.id[coord_num[i]] < t_min)
         {
            t_min = t.id[coord_num[i]];
            t_min -= 1.0*POW(10, -15);
            eqn_id = coord_num[i];
         }
      }

      cnd.node_id = 0;
      cnd.eqn_id = eqn_id;
      cnd.t_min = t_min;
      cnd.flag = 1;
      if(cnd.eqn_id == 0) cnd.N_OR_W = WEIGHT;
      else if(cnd.eqn_id > 0) cnd.N_OR_W = NODE;

      return cnd;
   }

}


static void ProjectNode(const Matrix eqn_matrix, const Vector dx, Vector x_projected)
{
   int i = -1, j = -1, l = -1, k = -1;
   int M = eqn_matrix.rows;
   int N = eqn_matrix.cols;
   double eqn_matrix_f[M*N];

   // convert to LAPACK memory layout
   MATRIX_TO_FORTRAN(M, N, eqn_matrix, eqn_matrix_f);

   // Perform QR
   int INFO = -1;
   int LDA = M;
   double TAU[MIN(M, N)];
   double WORK[N];
   dgeqr2_(&M, &N, eqn_matrix_f, &LDA, TAU, WORK, &INFO);

   double Q_EXPL[M][M], Q_TEMP[M][M];
   for(i = 0; i < M; ++i) for(j = 0; j < M; ++j) Q_TEMP[i][j] = 0.0;
   for(i = 0; i < M; ++i) Q_TEMP[i][i] = 1.0;


   // manually extract QR from LAPACK
   for(k = 0; k < N; ++k)
   {
      double H[M][M], v[M];
      for(i = 0; i < M; ++i) for(j = 0; j < M; ++j) Q_EXPL[i][j] = 0.0;
      for(i = 0; i < M; ++i) for(j = 0; j < M; ++j) H[i][j] = 0.0;


      for(i = 0; i < k; ++i) v[i] = 0.0;
      for(i = k+1; i < M; ++i) v[i] = eqn_matrix_f[i+k*M];
      v[k] = 1.0;


      for(i = 0; i < M; ++i)
         for(j = 0; j < M; ++j)
            H[i][j] = -TAU[k]*v[i]*v[j];

      for(i = 0; i < M; ++i) ++H[i][i];



      for(i = 0; i < M; ++i)
         for(j = 0; j < M; ++j)
            for(l = 0; l < M; ++l)
               Q_EXPL[i][j] += Q_TEMP[i][l]*H[l][j];


      for(i = 0; i < M; ++i)
         for(j = 0; j < M; ++j)
            Q_TEMP[i][j] = Q_EXPL[i][j];
   }


   /**********************************
   \* Obtain P = I - Q_RED x Q_RED^ \*
   **********************************/

   double Q_REDUCED[M][N];
   for(i = 0; i < M; ++i)
      for(j = 0; j < N; ++j)
         Q_REDUCED[i][j] = Q_EXPL[i][j];

   double Q_TRANS[N][M];
   for(i = 0; i < M; ++i) for(j = 0; j < N; ++j) Q_TRANS[j][i] = Q_REDUCED[i][j];


   Matrix PROJECTOR = Matrix_init(M, M);

   for(i = 0; i < M; ++i)
      for(j = 0; j < M; ++j)
         for(l = 0; l < N; ++l)
            PROJECTOR.id[i][j] += Q_REDUCED[i][l]*Q_TRANS[l][j];

   for(i = 0; i < M; ++i) for(j = 0; j < M; ++j) PROJECTOR.id[i][j] = -PROJECTOR.id[i][j];
   for(i = 0; i < M; ++i) ++PROJECTOR.id[i][i]; // add Identity


   memset(x_projected.id, 0, x_projected.len*sizeof(double));
   MatVec(PROJECTOR, dx, x_projected);

   Matrix_free(PROJECTOR);

   for(i = 0; i < M; ++i) x_projected.id[i] -= POW(10, -13);

}

#endif
