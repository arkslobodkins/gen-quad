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
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>

#define MAX_ELIM_WEIGHTS 10


static int CheckForFail(int info, double errorNorm, double errorNormPrev, Vector least_sq_sol);

static void ConstrVectDataReset(ConstrVectData *cVectData);
static ConstrNodeData ConstrNodeDataInit();
ATTR_UNUSED static void ConstrNodeDataReset(ConstrNodeData *cNodeData);

static int ConstrainedProjection(const quadrature *q_prev, quadrature *q_next);

static int ShortenVector(const RMatrix A, const Vector b, const quadrature *q_prev,
                         const quadrature *q_next, ConstrVectData *cVectData);

static ConstrNodeData ShortenNode(int node_num, const RMatrix A, const Vector b,
                                  const Vector z_old, const Vector z_new);

static int ProjectNode(const CMatrix eqn_matrix, const Vector dx, Vector x_projected);


// LeastSquaresNewton
// Receives initial quadrature guess. Primarily solves
// underdetermined systems of equations in the least
// squares sense. Returns success if algorithm converged
// and all nodes are inside of the domain and if all nodes are positive.
bool LeastSquaresNewton(const bool_enum FLAG_CONSTR, const int_fast8_t *basis, quadrature *q_orig, int *its)
{
   assert(q_orig->k >= 0);

   int elim_weights = 0;
   ConstrVectData cVectData = ConstrVectDataInit();

   int k               = q_orig->k;
   int deg             = q_orig->deg;
   int numFuncs        = q_orig->num_funcs;
   int dim             = q_orig->dim;
   int *dims           = q_orig->dims;
   bool SOL_FLAG       = SOL_NOT_FOUND;

   int itsLoc       = 0; int maxiter = 25;
   double errorNorm = -1.0, errorNormPrev = -1.0, errorNormUpdate = -1;
   double q_tol     = QUAD_TOL; // 10^(-15);

   // initialize LAPACK
   int nrows             = numFuncs, ncols = (dim+1)*k;
   CMatrix JACOBIAN      = CMatrix_init(nrows, ncols);
   Vector RHS            = Vector_init(numFuncs);
   char TRANS = 'N';
   int INFO = 0, NRHS = 1, LDA = nrows, LEAD_DIM = MAX(nrows, ncols), LWORK = 5*ncols, SMALL_DIM = MIN(nrows, ncols);
   double *WORK        = (double *)malloc(LWORK*size_double);
   Vector LEAST_SQ_SOL = Vector_init(LEAD_DIM);

   quadrature *q_prev = quadrature_make_full_copy(q_orig);
   quadrature *q_next = quadrature_make_full_copy(q_orig);

   // return if input is a satisfactory quadrature
   GetFunction(basis, q_prev, RHS);
   errorNorm = V_ScaledTwoNorm(RHS);
   if( (errorNorm < q_tol) && (QuadInConstraint(q_prev) == true) )
   {
      SOL_FLAG = SOL_FOUND;
      *its = 0;
      goto FREERETURN;
   }


   do
   {

      if(FLAG_CONSTR == ON)
      {
         if(cVectData.ACTIVE == ON)
         {
            if(cVectData.N_OR_W == WEIGHT)
            {
               if(k == 1  || ++elim_weights > MAX_ELIM_WEIGHTS)
               {
                  SOL_FLAG = SOL_NOT_FOUND;
                  goto FREERETURN;
               }
               ncols = (dim+1)*--k;
               SMALL_DIM = MIN(nrows, ncols);
               LEAD_DIM  = MAX(nrows, ncols);
               quadrature_remove_element(cVectData.boundaryNodeId, q_next);
               quadrature_remove_element(cVectData.boundaryNodeId, q_prev);
            }
         }
      }

      GetJacobian(basis, q_prev, JACOBIAN);
      GetFunction(basis, q_prev, RHS);

      for(int i = 0; i < SMALL_DIM; ++i)        LEAST_SQ_SOL.id[i] = RHS.id[i];
      for(int i = SMALL_DIM; i < LEAD_DIM; ++i) LEAST_SQ_SOL.id[i] = 0.0;

      // solve least squares problem using LAPACK routine and time it
      time_t start = clock();
      dgels_(&TRANS, &nrows, &ncols, &NRHS, JACOBIAN.id, &LDA,
             &LEAST_SQ_SOL.id[0], &LEAD_DIM, WORK, &LWORK, &INFO);
      time_t end = clock(); LSQ_TIME += (double)(end - start)/CLOCKS_PER_SEC;

      for(int i = 0; i < ncols; ++i) q_next->z.id[i] = q_prev->z.id[i] - LEAST_SQ_SOL.id[i];

      errorNormPrev = errorNorm;
      GetFunction(basis, q_next, RHS);
      errorNorm = V_ScaledTwoNorm(RHS);
      int check_values = CheckForFail(INFO, errorNorm, errorNormPrev, LEAST_SQ_SOL);
      if(check_values != Q_SUCCESS)
      {
         SOL_FLAG = SOL_NOT_FOUND;
         goto FREERETURN;
      }

      if(FLAG_CONSTR == ON)
      {
         int P_FLAG = ConstrainedProjection(q_prev, q_next);
         if(P_FLAG != CONSTR_SUCCESS)
         {
            SOL_FLAG = SOL_NOT_FOUND;
            goto FREERETURN;
         }
         int C_FLAG = ConstrainedOptimization(q_prev, q_next, &cVectData);
         if(C_FLAG != CONSTR_SUCCESS)
         {
            SOL_FLAG = SOL_NOT_FOUND;
            goto FREERETURN;
         }
      }
      else
      {
         if( !QuadInConstraint(q_next) && itsLoc > 5)
         {
            SOL_FLAG = SOL_NOT_FOUND;
            goto FREERETURN;
         }
      }

      GetFunction(basis, q_next, RHS);
      errorNormUpdate = V_ScaledTwoNorm(RHS);

      quadrature_assign(q_next, q_prev);
      ++itsLoc;

   } while( (itsLoc < maxiter) && (errorNormUpdate > q_tol) );


   // check if quadrature satisfies constraints
   if( ( !QuadInConstraint( q_next ))
    || (  errorNormUpdate > q_tol) )
   {
      SOL_FLAG = SOL_NOT_FOUND;
      goto FREERETURN;
   }

   // check whether Newton's method succeeded
   if( (itsLoc < maxiter)
    && (isnan(errorNormUpdate) == 0)
    && (isinf(errorNormUpdate) == 0)
    && (errorNormUpdate <= q_tol) )
   {
      quadrature_realloc(q_next->k, dim, dims, deg, q_orig);
      quadrature_assign(q_next, q_orig);
      SOL_FLAG = SOL_FOUND;
   }
   else
      SOL_FLAG = SOL_NOT_FOUND;


FREERETURN:
   free(WORK);
   Vector_free(LEAST_SQ_SOL);
   Vector_free(RHS);
   CMatrix_free(JACOBIAN);
   quadrature_free(q_prev);
   quadrature_free(q_next);

   *its = itsLoc;
   return SOL_FLAG;
}// end LeastSquaresNewton



static int CheckForFail(int INFO, double errorNorm, double errorNormPrev, Vector least_sq_sol)
{

   if(errorNorm > errorNormPrev+2)     // fail if method is not converging
      return NOT_CONVERGE;

   if(V_InfNorm(least_sq_sol) > 100)   // fail if solution is exploding
      return DIVERGE_ERR;

   else if(INFO != 0)                  // fail if LAPACK routine has failed
      return LAPACK_ERR;

   else if(V_CheckInf(least_sq_sol))
      return INF_VAL;

   else if(V_CheckNan(least_sq_sol))
      return NAN_VAL;

   else return Q_SUCCESS;
}



ConstrVectData ConstrVectDataInit()
{
   ConstrVectData cVectData;
   cVectData.boundaryNodeId = -1;
   cVectData.eqnId          = -1;
   cVectData.tMin           =  1.0;
   cVectData.ACTIVE         = OFF;
   cVectData.N_OR_W         = NONE;

   return cVectData;
}

static void ConstrVectDataReset(ConstrVectData *cVectData)
{
   cVectData->boundaryNodeId = -1;
   cVectData->eqnId          = -1;
   cVectData->tMin           =  1.0;
   cVectData->ACTIVE         = OFF;
   cVectData->N_OR_W         = NONE;
}

static ConstrNodeData ConstrNodeDataInit()
{
   ConstrNodeData cNodeData;
   cNodeData.nodeId = -1;
   cNodeData.eqnId  = -1;
   cNodeData.tMin   =  1.0;
   cNodeData.ACTIVE =  OFF;
   cNodeData.N_OR_W =  NONE;

   return cNodeData;
}

ATTR_UNUSED static void ConstrNodeDataReset(ConstrNodeData *cNodeData)
{
   cNodeData->nodeId = -1;
   cNodeData->eqnId  = -1;
   cNodeData->tMin   =  1.0;
   cNodeData->ACTIVE =  OFF;
   cNodeData->N_OR_W =  NONE;
}


static int ConstrainedProjection(const quadrature *q_prev, quadrature *q_next)
{
   int i,j,d;
   int k = q_prev->k;
   int dim = q_prev->dim;
   int RET_FLAG = CONSTR_SUCCESS;
   assert(dim == q_next->constr->M.cols);

   int num_eqns = q_next->constr->M.rows;
   RMatrix M = q_next->constr->M;

   // allocate vectors on the stack
   Vector node_change; node_change.len = dim; double node_change_id[dim]; node_change.id = node_change_id;
   Vector node_projected; node_projected.len = dim; double node_projected_id[dim]; node_projected.id = node_projected_id;
   memset(node_change.id, 0, dim*sizeof(double)); memset(node_projected.id, 0, dim*sizeof(double));

   for(i = 0; i < k; ++i)
   {
      int active_eqn_count = 0, node_index = i*dim;
      bool do_project = false;
      bool eqn_flags[num_eqns]; for(j = 0; j < num_eqns; ++j) eqn_flags[j] = false;

      for(j = 0; j < num_eqns; ++j)
      {
         if(QuadEqnOnTheBoundary(q_prev, i, j))
         {
            ++active_eqn_count;
            eqn_flags[j] = true;
            if( !do_project)
               if( !QuadInDomainElem(q_next, i) ) do_project = true;
         }
      }

      if(active_eqn_count > 0  && do_project == true)
      {
         int count = 0;
         CMatrix eqn_matrix = CMatrix_init(dim, active_eqn_count);

         for(j = 0; j < num_eqns; ++j)
         {
            if(eqn_flags[j] == true)
            {
               for(d = 0; d < dim; ++d)
                  eqn_matrix.cid[count][d]= M.rid[j][d];

               ++count;
            }
         }

         for(d = 0; d < dim; ++d)
            node_change.id[d] = q_next->x[node_index+d] - q_prev->x[node_index+d];

         int P_FLAG = ProjectNode(eqn_matrix, node_change, node_projected);
         if(P_FLAG != CONSTR_SUCCESS)
         {
            CMatrix_free(eqn_matrix);
            RET_FLAG = CONSTR_FAILURE;
            return RET_FLAG;
         }
         for(d = 0; d < dim; ++d)
            q_next->x[node_index+d] = q_prev->x[node_index+d] + node_projected.id[d] - POW_DOUBLE(10, -13);

         CMatrix_free(eqn_matrix);
      }

   }
   COND_TEST_1;

   return RET_FLAG;
}


int ConstrainedOptimization(const quadrature *q_prev, quadrature *q_next, ConstrVectData *cVectData)
{
   assert(q_prev->k == q_next->k);

   int q_len =  q_next->z.len;
   int RET_FLAG = CONSTR_SUCCESS;

   Vector q_diff = Vector_init(q_len);
   quadrature *q_next_copy = quadrature_make_full_copy(q_next);
   for(int i = 0; i < q_len; ++i)
      q_diff.id[i] = q_prev->z.id[i] - q_next_copy->z.id[i];

   if( ( QuadInConstraint( q_next_copy) == false )
         && ( QuadInConstraint(q_prev ) == true ) )
   {
      RET_FLAG = ShortenVector(q_prev->constr->M_FULL, q_prev->constr->b_FULL,
                               q_prev, q_next_copy, cVectData);
      COND_TEST_2;

      // return if quadrature vector was not shortened successfully
      if(RET_FLAG != CONSTR_SUCCESS && cVectData->ACTIVE == OFF)
      {
         ConstrVectDataReset(cVectData);
         Vector_free(q_diff);
         return RET_FLAG;
      }

      for(int i = 0; i < q_len; ++i)
         q_next_copy->z.id[i] = q_prev->z.id[i] + (-cVectData->tMin + POW_DOUBLE(10.0, -14) ) * q_diff.id[i];

      COND_TEST_3;

      if(QuadInConstraint( q_next_copy ) == 1)
      {
         quadrature_assign(q_next_copy, q_next);
         RET_FLAG = CONSTR_SUCCESS;
      }
      else
         RET_FLAG = CONSTR_FAILURE;

   }
   else
   {
      ConstrVectDataReset(cVectData);
      RET_FLAG = CONSTR_SUCCESS;
   }

   quadrature_free(q_next_copy);
   Vector_free(q_diff);
   return RET_FLAG;
}

// Shortens q_next vector such that every node satisfies inequality  A*q_next_i = b,
// provided that q_prev satisfies the constraints.
static int ShortenVector(const RMatrix A, const Vector b, const quadrature *q_prev, const quadrature *q_next, ConstrVectData *cVectData)
{
   assert(A.cols == q_prev->dim+1);
   assert(q_prev->k == q_next->k);

   int i, count;
   int k = q_prev->k;
   int dim = q_prev->dim;
   int *eqn = (int *)calloc(k, sizeof(int));

   ConstrNodeData *cNodeDataArr = (ConstrNodeData *)malloc(k* sizeof(ConstrNodeData));
   for(int i = 0; i < k; ++i)
      cNodeDataArr[i] = ConstrNodeDataInit();

   // allocate vectors on the stack
   Vector z_prev_node = {0}; z_prev_node.len = dim+1;
   double z_prev_id[dim+1]; z_prev_node.id = z_prev_id;
   Vector z_next_node = {0}; z_next_node.len = dim+1;
   double z_next_id[dim+1]; z_next_node.id = z_next_id;

   for(count = 0, i = 0; i < k; ++i)
   {
      quadrature_get_elem(q_prev, i, z_prev_node);
      quadrature_get_elem(q_next, i, z_next_node);

      if( !QuadInConstraintElem(q_next, i) )
      {
         cNodeDataArr[i] = ShortenNode(i, A, b, z_prev_node, z_next_node);
         if(cNodeDataArr[i].ACTIVE == ON)
            ++count;
         else if(cNodeDataArr[i].ACTIVE == OFF)
         {
            ConstrVectDataReset(cVectData);
            free(cNodeDataArr);
            free(eqn);
            return CONSTR_FAILURE;
         }
      }
   }

   Vector tVec = Vector_init(k);
   for(int i = 0; i < k; ++i)
      tVec.id[i] = cNodeDataArr[i].tMin;

   if(count > 0)
   {
      VMin vMin = VectorMin(tVec);
      int min_index             = vMin.min_index;
      cVectData->boundaryNodeId = min_index;
      cVectData->eqnId          = cNodeDataArr[min_index].eqnId;
      cVectData->tMin           = vMin.min_value;
      cVectData->ACTIVE         = ON;
      cVectData->N_OR_W         = cNodeDataArr[min_index].N_OR_W;
   }
   else
      ConstrVectDataReset(cVectData);

   Vector_free(tVec);
   free(cNodeDataArr);
   free(eqn);

   return CONSTR_SUCCESS;
}


// Returns data information needed to map z_new onto the boundary, such that A*z_new = b_bound,
// provided that z_old satisfies A*z_old <= b_bound, and A*z_new > b_bound.
// If both A*z_new <= b_bound, and A*z_old <= b_bound,
// the routine does no further computations and default is returned.
static ConstrNodeData ShortenNode(int node_num, const RMatrix A, const Vector b_bound, const Vector z_old, const Vector z_new)
{
   assert(A.rows == b_bound.len);

   int i;
   int out_count;
   int nrows = A.rows;
   int ncols = A.cols;
   ConstrNodeData cNodeData = ConstrNodeDataInit();

   // Allocate vectors on the stack
   Vector b_old, b_new, b_diff, t, dz;
   b_old.len = nrows; b_new.len = nrows; b_diff.len = nrows; t.len = nrows; dz.len = ncols;
   double b_old_id[nrows], b_new_id[nrows], b_diff_id[nrows], t_id[nrows], dz_id[ncols];
   b_old.id = b_old_id; b_new.id = b_new_id; b_diff.id = b_diff_id; t.id = t_id; dz.id = dz_id;
   memset(t.id, 0, nrows*sizeof(double));
   memset(b_old.id, 0, nrows*sizeof(double));
   memset(b_diff.id, 0, nrows*sizeof(double));
   memset(b_new.id, 0, nrows*sizeof(double));
   int outEqn[nrows];
   memset(outEqn, -1, nrows*sizeof(int));


   VectorAddScale(1.0, z_new, -1.0, z_old, dz);
   RMatVec(A, dz, b_diff);
   RMatVec(A, z_new, b_new);
   RMatVec(A, z_old, b_old);
   // return if z_old does not satisfy inequality constraints
   for (i = 0; i < nrows; ++i)
      if(b_old.id[i] > b_bound.id[i])
         return cNodeData;

   // compute t for those equations that z_new does not satisfy inequality constraints
   for (out_count = 0, i = 0; i < nrows; ++i)
      if (b_new.id[i] > b_bound.id[i])
      {
         t.id[out_count] = (b_bound.id[i] - b_old.id[i]) / b_diff.id[i];
         outEqn[out_count++] = i;
      }

   // return if z_new satisfies constraints
   if(out_count == 0)
      return cNodeData;
   // compute and return data if z_new does not satisfy constraints
   else
   {
      int eqnId = outEqn[0];
      double tMin = t.id[0];
      for(i = 1; i < out_count; ++i)
         if(t.id[i] < tMin)
         {
            tMin  = t.id[i];
            eqnId = outEqn[i];
         }

      cNodeData.nodeId = node_num;
      cNodeData.eqnId  = eqnId;
      cNodeData.tMin   = tMin;
      cNodeData.ACTIVE = ON;
      if(cNodeData.eqnId == 0)     cNodeData.N_OR_W = WEIGHT;
      else if(cNodeData.eqnId > 0) cNodeData.N_OR_W = NODE;

      return cNodeData;
   }

}


// Projects dx onto equations specified by eqn_matrix. Reduced matrix Q of eqn_matrix
// is extracted from QR factorization, and projector P is computed by
// P = I-Q_red*Q_red'. Results are stored in x_projected.
static int ProjectNode(const CMatrix eqn_matrix, const Vector dx, Vector x_projected)
{
   int i, j, l, k;

   // convert to LAPACK memory layout
   int M = eqn_matrix.rows;
   int N = eqn_matrix.cols;

   // Perform QR
   int INFO = -1, LDA = M;
   double TAU[MIN(M, N)], WORK[N];
   dgeqr2_(&M, &N, eqn_matrix.id, &LDA, TAU, WORK, &INFO);

   if(INFO != 0)
      return LAPACK_ERR;

   double Q_EXPL[M][M], Q_TEMP[M][M];
   for(i = 0; i < M; ++i) for(j = 0; j < M; ++j) Q_TEMP[i][j] = 0.0;
   for(i = 0; i < M; ++i) Q_TEMP[i][i] = 1.0;

   // manually extract QR from LAPACK
   for(k = 0; k < N; ++k)
   {
      double H[M][M], v[M];
      for(i = 0; i < M; ++i) for(j = 0; j < M; ++j) Q_EXPL[i][j] = 0.0;
      for(i = 0; i < M; ++i) for(j = 0; j < M; ++j) H[i][j] = 0.0;

      for(i = 0; i < k; ++i)   v[i] = 0.0;
      for(i = k+1; i < M; ++i) v[i] = C_ELEM_ID(eqn_matrix, i, k);
      v[k] = 1.0;

      for(i = 0; i < M; ++i)
         for(j = 0; j < M; ++j)
            H[i][j] = -TAU[k]*v[i]*v[j];

      for(i = 0; i < M; ++i) ++H[i][i];

      for(i = 0; i < M; ++i)
         for(l = 0; l < M; ++l)
            for(j = 0; j < M; ++j)
               Q_EXPL[i][j] += Q_TEMP[i][l]*H[l][j];

      for(i = 0; i < M; ++i)
         for(j = 0; j < M; ++j)
            Q_TEMP[i][j] = Q_EXPL[i][j];
   }


   /******************************************
   \* Obtain P = I - Q_REDUCED x Q_REDUCED^ \*
   ******************************************/
   double Q_REDUCED[M][N];
   for(i = 0; i < M; ++i)
      for(j = 0; j < N; ++j)
         Q_REDUCED[i][j] = Q_EXPL[i][j];

   double Q_TRANS[N][M];
   for(i = 0; i < M; ++i)
      for(j = 0; j < N; ++j)
         Q_TRANS[j][i] = Q_REDUCED[i][j];

   RMatrix PROJECTOR = RMatrix_init(M, M);
   for(i = 0; i < M; ++i)
      for(l = 0; l < N; ++l)
         for(j = 0; j < M; ++j)
            PROJECTOR.rid[i][j] += Q_REDUCED[i][l]*Q_TRANS[l][j];

   for(i = 0; i < M; ++i)
      for(j = 0; j < M; ++j)
         PROJECTOR.rid[i][j] = -PROJECTOR.rid[i][j];
   for(i = 0; i < M; ++i) ++PROJECTOR.rid[i][i];

   memset(x_projected.id, 0, x_projected.len*size_double);
   RMatVec(PROJECTOR, dx, x_projected);

   RMatrix_free(PROJECTOR);

   return CONSTR_SUCCESS;
}

