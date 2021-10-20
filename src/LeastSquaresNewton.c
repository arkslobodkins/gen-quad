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


static bool CheckForFail(int info, double error_norm, double error_norm_prev, Vector least_sq_sol);

#ifdef CONSTR_OPT
static ConstrVectData ConstrVectDataInit();
static void ConstrVectDataReset(ConstrVectData  *cVectData);
static ConstrNodeData ConstrNodeDataInit();
static void ConstrNodeDataReset(ConstrNodeData  *cNodeData);

static void ConstrainedProjection(const_quadrature *q_prev, quadrature *q_next);
static bool ConstrainedOptimization(const_quadrature *q_prev, quadrature *q_next, ConstrVectData *cVecData);
static ConstrNodeData ShortenNode(const Matrix A, const Vector b, const Vector z_old, const Vector z_new);
static void ProjectNode(const Matrix eqn_matrix, const Vector dx, Vector x_projected);
#endif




/* LeastSquaresNewton
 * Receives initial quadrature guess. Primarily solves
 * underdetermined systems of equations in the least
 * squares sense. Returns success if algorithm converged
 * and all nodes are inside of the domain and if all nodes are positive..
 */
bool LeastSquaresNewton(const bool_enum FLAG_CONSTR, const int_fast8_t *basis, quadrature *q_orig, int *its)
{
   assert(q_orig->k >= 0);
   int k = q_orig->k;

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
   quad_set_funcs_and_constr(q_prev);

   quadrature *q_next = quadrature_init(k, dim, dims, p, D);
   quad_set_funcs_and_constr(q_next);

   quadrature_assign(q_orig, q_prev);
   quadrature_assign(q_orig, q_next);

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
   ConstrVectData cVectData = ConstrVectDataInit();
#else
   if(FLAG_CONSTR == ON)
      PRINT_ERR("constrained optimization is turned off and won't be performed", __LINE__, __FILE__);

#endif


   // return if input is a satisfactory quadrature
   GetFunction(basis, *q_prev, rhs.id);
   error_norm = V_ScaledTwoNorm(rhs);
   if( (error_norm < tol) && (q_prev->inConstraint( (const_quadrature *)q_prev ) == true) )
   {
      SOL_FLAG = SOL_FOUND;
      *its = 0;
      goto FREERETURN;
   }


   do
   {

#ifdef CONSTR_OPT
      if(cVectData.ACTIVE_CONSTRAINTS == ON)
      {
         switch(cVectData.N_OR_W)
         {
            case WEIGHT:
               ++elim_weights;
               if(k == 1  || elim_weights > MAX_ELIM_WEIGHTS)
               {
                  SOL_FLAG = SOL_NOT_FOUND;
                  goto FREERETURN;
               }

               --k;
               n_cols = (dim+1)*k;
               n_rows = num_funcs;
               quadrature_remove_element(cVectData.boundaryNodeId, q_next);
               quadrature_remove_element(cVectData.boundaryNodeId, q_prev);
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
      GetJacobian(basis, *q_prev, jacobian.id);
      GetFunction(basis, *q_prev, rhs.id);

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
      GetFunction(basis, *q_next, rhs.id);
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
         ConstrainedProjection((const_quadrature *)q_prev, q_next);
         CONSTR_RETURN = ConstrainedOptimization((const_quadrature *)q_prev, q_next, &cVectData);

         if(CONSTR_RETURN == CONSTRAINT_FAILURE)
         {
            SOL_FLAG = SOL_NOT_FOUND;
            goto FREERETURN;
         }
      }
#endif


      GetFunction(basis, *q_next, rhs.id);
      error_norm_update = V_ScaledTwoNorm(rhs);
      quadrature_assign(q_next, q_prev);
      ++its_loc;

   } while( (its_loc < maxiter) && (error_norm_update > tol) );


   // check if quadrature satisfies constraints
   if( (q_next->inConstraint( (const_quadrature *)q_next ) == false) || (error_norm_update > tol) )
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
      quadrature_assign(q_next, q_orig);
      *its = its_loc;
      SOL_FLAG = SOL_FOUND;
   }
   else
   {
      *its = its_loc;
      SOL_FLAG = SOL_NOT_FOUND;
   }

FREERETURN:
   Vector_free(least_sq_sol);
   Vector_free(work);
   Vector_free(rhs);
   Vector_free(jacobian_transpose);
   Vector_free(jacobian);
   quadrature_free(q_prev);
   quadrature_free(q_next);

   return SOL_FLAG;
}// end LeastSquaresNewton



static bool CheckForFail(int info, double error_norm, double error_norm_prev, Vector least_sq_sol)
{

   if( error_norm > (error_norm_prev+10) )      // return if Newton's method does not have reasonable convergence
      return 1;

   else if( V_InfNorm(least_sq_sol) > 100 )    // return if solution is exploding
      return 1;

   else if(info > 0)                           // return if LAPACK routine has failed
      return 1;

   else if( V_CheckInfAndNan(least_sq_sol) )
      return 1;

   else return 0;

}


#ifdef CONSTR_OPT

static ConstrVectData ConstrVectDataInit()
{
   ConstrVectData cVectData;
   cVectData.ACTIVE_CONSTRAINTS = OFF;
   cVectData.N_OR_W             = NONE;
   cVectData.boundaryNodeId     = -1;
   cVectData.eqnId              = -1;

   return cVectData;
}

static void ConstrVectDataReset(ConstrVectData *cVectData)
{
   cVectData->ACTIVE_CONSTRAINTS = OFF;
   cVectData->N_OR_W             = NONE;
   cVectData->boundaryNodeId     = -1;
   cVectData->eqnId              = -1;

}


static ConstrNodeData ConstrNodeDataInit()
{
   ConstrNodeData cNodeData;
   cNodeData.eqnId  = -1;
   cNodeData.nodeId = -1;
   cNodeData.tMin   = 1.0;
   cNodeData.flag    = 0;
   cNodeData.N_OR_W  = NONE;

   return cNodeData;
}

static void ConstrNodeDataReset(ConstrNodeData *cNodeData)
{
   cNodeData->eqnId  = -1;
   cNodeData->nodeId = -1;
   cNodeData->tMin   = 1.0;
   cNodeData->flag   = 0;
   cNodeData->N_OR_W = NONE;
}


static void ConstrainedProjection(const_quadrature *q_prev, quadrature *q_next)

{
   int k = q_prev->k;
   int dim = q_prev->params->dim;
   assert(dim == q_next->cons->M.cols);
   double tol = pow(10, -12);

   int num_eqns = q_next->cons->M.rows;
   Matrix M = q_next->cons->M;
   Vector b = q_next->cons->b;

   Vector node_change = Vector_init(dim);
   Vector node_projected = Vector_init(dim);

   for(int i = 0; i < k; ++i)
   {
      int active_eqn_count = 0, node_index = i*dim;
      bool do_project = false;
      bool eqn_flags[num_eqns]; for(int n = 0; n < num_eqns; ++n) eqn_flags[n] = false;

      for(int j = 0; j < num_eqns; ++j)
      {
         double lhs_prev = 0.0;
         for(int d = 0; d < dim; ++d)
            lhs_prev += M.id[j][d] * q_prev->x[node_index+d];

         if( fabs( lhs_prev - b.id[j] ) <= tol )
         {
            ++active_eqn_count;
            eqn_flags[j] = true;
            if( ! do_project)
               if( ! q_next->inDomainElem((const_quadrature *)q_next, i) ) do_project = true;
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


bool ConstrainedOptimization(const_quadrature *q_prev, quadrature *q_next, ConstrVectData *cVectData)
{
   assert(q_prev->k == q_next->k);

   int k = q_prev->k;
   int dim = q_prev->params->dim;
   int q_len = (dim+1)*k;

   Vector q_diff = Vector_init(q_len);
   for(int i = 0; i < q_len; ++i)
      q_diff.id[i] = q_prev->z[i] - q_next->z[i];

   // Compute only if equations at the previous iteration satisfy constraints
   if( ( q_next->inConstraint( (const_quadrature *)q_next) == false )
         && ( q_prev->inConstraint(q_prev ) == true ) )
   {
      ConstrNodeData cNodeData = constrain_vector(q_prev->FULL_A, q_prev->FULL_b,
                                                  q_prev, (const_quadrature *)q_next);
      COND_TEST_2;

      if(cNodeData.flag == 0) {
         ConstrVectDataReset(cVectData);
         Vector_free(q_diff);
         return CONSTRAINT_FAILURE;
      }

      switch(cNodeData.N_OR_W)
      {
      case NODE:
         for(int i = 0; i < q_len; ++i)
            q_next->z[i] = q_prev->z[i] + (-cNodeData.tMin + POW(10, -12) ) * q_diff.id[i];

         cVectData->ACTIVE_CONSTRAINTS = ON;
         cVectData->N_OR_W = NODE;
         cVectData->boundaryNodeId = cNodeData.nodeId;
         cVectData->eqnId = cNodeData.eqnId;
         break;

      case WEIGHT:
         for(int i = 0; i < q_len; ++i)
            q_next->z[i] = q_prev->z[i] + (-cNodeData.tMin + POW(10, -12) ) * q_diff.id[i];

         cVectData->ACTIVE_CONSTRAINTS = ON;
         cVectData->N_OR_W = WEIGHT;
         cVectData->boundaryNodeId = cNodeData.nodeId;
         cVectData->eqnId = cNodeData.eqnId;
         break;

      default:
         ConstrVectDataReset(cVectData);
         break;
      }
      COND_TEST_3;

      if(q_next->inConstraint( (const_quadrature *)q_next) == 1) {
         Vector_free(q_diff);
         return CONSTRAINT_SUCCESS;
      }
      else {
         Vector_free(q_diff);
         return CONSTRAINT_FAILURE;
      }
   }
   else if( ( q_next->inConstraint( (const_quadrature *)q_next) == true )
              && ( q_prev->inConstraint(q_prev ) == true ) ) {
      ConstrVectDataReset(cVectData);
      Vector_free(q_diff);
      return CONSTRAINT_SUCCESS;
   }
   else {
      ConstrVectDataReset(cVectData);
      Vector_free(q_diff);
      return CONSTRAINT_FAILURE;
   }

}


ConstrNodeData constrain_vector(const Matrix A, const Vector b, const_quadrature *q_prev, const_quadrature *q_next)
{
   assert(A.cols == q_prev->params->dim+1);
   assert(q_prev->k == q_next->k);

   int k = q_prev->k;
   int dim = q_prev->params->dim;
   ConstrNodeData *cNodeDataArr = (ConstrNodeData *)malloc( k* sizeof(ConstrNodeData) );

   int *eqn = (int *)calloc( k, sizeof(int) );
   for(int i = 0; i < k; ++i)
      cNodeDataArr[i] = ConstrNodeDataInit();

   // allocate vectors on the stack
   Vector z_prev_node = {0}; z_prev_node.len = dim+1; double z_prev_id[dim+1]; z_prev_node.id = z_prev_id;
   Vector z_next_node = {0}; z_next_node.len = dim+1; double z_next_id[dim+1]; z_next_node.id = z_next_id;

   int i, count;
   ConstrNodeData cNodeData;
   for(count = 0, i = 0; i < k; ++i) {
      quadrature_get_elem(q_prev, i, z_prev_node);
      quadrature_get_elem(q_next, i, z_next_node);

      if(q_next->inConstraintElem(q_next, i) == false) {
         cNodeDataArr[i] = ShortenNode(A, b, z_prev_node, z_next_node);
         if(cNodeDataArr[i].flag == 1) // store equations if succeeded
            ++count;
         else if(cNodeDataArr[i].flag == 0)
         {
            ConstrNodeDataReset(&cNodeData);
            goto FREERETURN;
         }
      }
   }

   Vector tVec = Vector_init(k);
   for(int i = 0; i < k; ++i)
      tVec.id[i] = cNodeDataArr[i].tMin;

   int min_index = -1;
   VMin vMin = {0};
   if(count > 0) // compute minimum t for mapping nodes to the boundary, defaults to 1
   {
      vMin = VectorMin(tVec);
      min_index        = vMin.min_index;
      cNodeData.N_OR_W = cNodeDataArr[min_index].N_OR_W;
      cNodeData.eqnId  = cNodeDataArr[min_index].eqnId;
      cNodeData.nodeId = min_index;
      cNodeData.tMin   = vMin.min_value;
      cNodeData.flag   = 1;
   }
   else
      ConstrNodeDataReset(&cNodeData);

FREERETURN:
   Vector_free(tVec);
   free(cNodeDataArr);
   free(eqn);

   return cNodeData;
}


static ConstrNodeData ShortenNode(const Matrix A, const Vector b_bound, const Vector z_old, const Vector z_new)
{
   assert(A.rows == b_bound.len);

   int i = -1, j = -1;
   int out_count = -1;
   int nrows = A.rows;
   int ncols = A.cols;
   int bound_eqn[nrows];
   ConstrNodeData cNodeData;

   // Allocate vectors on the stack
   Vector b_old, b_new, b_diff, t, dz;
   b_old.len = nrows; b_new.len = nrows; b_diff.len = nrows; t.len = nrows; dz.len = ncols;
   double b_old_id[nrows], b_new_id[nrows], b_diff_id[nrows], t_id[nrows], dz_id[ncols];
   b_old.id = b_old_id; b_new.id = b_new_id; b_diff.id = b_diff_id; t.id = t_id; dz.id = dz_id;
   memset(t.id, 0, nrows*sizeof(double));
   memset(b_old.id, 0, nrows*sizeof(double));
   memset(b_diff.id, 0, nrows*sizeof(double));
   memset(b_new.id, 0, nrows*sizeof(double));
   memset(bound_eqn, -1, nrows*sizeof(int));


   VectorAddScale(1.0, z_new, -1.0, z_old, dz);
   MatVec(A, dz, b_diff);
   MatVec(A, z_new, b_new);
   MatVec(A, z_old, b_old);
   for (i = 0; i < nrows; ++i) {
      if(b_old.id[i] > b_bound.id[i]) {
         ConstrNodeDataReset(&cNodeData);
         return cNodeData;
      }
   }

   for (out_count = 0, i = 0; i < nrows; ++i) {
      if ( b_new.id[i] > b_bound.id[i] ) {
         t.id[i] = (b_bound.id[i] - b_old.id[i]) / b_diff.id[i];
         bound_eqn[out_count] = i;
         ++out_count;
      }
   }


   if(out_count == 0)  // exit if all nodes already satisfy constraints
   {
      ConstrNodeDataReset(&cNodeData);
      return cNodeData;
   }
   else // compute minimum t
   {

      int eqnId = bound_eqn[0];
      double tMin = t.id[bound_eqn[0]];
      for(i = 1; i < out_count; ++i) {
         if (t.id[bound_eqn[i]] < tMin) {
            tMin = t.id[bound_eqn[i]];
            eqnId = bound_eqn[i];
         }
      }
      cNodeData.nodeId = 0;
      cNodeData.eqnId  = eqnId;
      cNodeData.tMin   = tMin;
      cNodeData.flag    = 1;
      if(cNodeData.eqnId == 0) cNodeData.N_OR_W = WEIGHT;
      else if(cNodeData.eqnId > 0) cNodeData.N_OR_W = NODE;

      return cNodeData;
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
