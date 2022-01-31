#include "ConstrainedOptimization.h"
#include "Quadrature.h"
#include "Constraints.h"
#include "Conditional_Debug.h"
#include "LINALG.h"

#include <assert.h>
#include <string.h>
#include <stdlib.h>


int ConstrainedProjection(const quadrature *q_prev, quadrature *q_next)
{
   int i,j,d;
   int k = q_prev->num_nodes;
   int dim = q_prev->dim;
   assert(dim == q_next->constr->M.cols);

   int num_eqns = q_next->constr->M.rows;
   RMatrix M = q_next->constr->M;

   StaticVectorInit(node_change, dim);
   StaticVectorInit(node_projected, dim);

   for(i = 0; i < k; ++i)
   {
      int active_eqn_count = 0;
      int node_index = i*dim;
      bool do_project = false;
      bool eqn_flags[num_eqns];
      for(j = 0; j < num_eqns; ++j) eqn_flags[j] = false;

      if( QuadOnTheBoundary(q_prev, i)  && !QuadInDomainElem(q_next, i) )
         do_project = true;

      if(do_project)
      {
         for(j = 0; j < num_eqns; ++j)
         {
            if(QuadEqnOnTheBoundary(q_prev, i, j)) {
               ++active_eqn_count;
               eqn_flags[j] = true;
            }
         }
         CMatrix eqn_matrix = CMatrix_init(dim, active_eqn_count);
         int count = 0;
         for(j = 0; j < num_eqns; ++j) {
            if(eqn_flags[j] == true) {
               for(d = 0; d < dim; ++d)
                  eqn_matrix.cid[count][d]= M.rid[j][d];
               ++count;
            }
         }
         for(d = 0; d < dim; ++d)
            node_change.id[d] = q_next->x[node_index+d] - q_prev->x[node_index+d];
         int P_FLAG = ProjectNode(eqn_matrix, node_change, node_projected);
         if(P_FLAG != CONSTR_SUCCESS) {
            CMatrix_free(eqn_matrix);
            return P_FLAG;
         }
         for(d = 0; d < dim; ++d)
            q_next->x[node_index+d] = q_prev->x[node_index+d] + node_projected.id[d] - POW_DOUBLE(10, -13);

         CMatrix_free(eqn_matrix);
      }
   }
   COND_TEST_1;
   return CONSTR_SUCCESS;
}


int ConstrainedOptimization(ConstrOptData *data, const quadrature *q_prev, quadrature *q_next, ConstrVectData *cVectData)
{
   assert(q_prev->num_nodes == q_next->num_nodes);
   int qlen = q_next->z.len;

   Vector q_diff           = data->q_diff;
   quadrature *q_next_copy = data->q_next_copy;
   for(int i = 0; i < qlen; ++i)
      q_diff.id[i] = q_next_copy->z.id[i] - q_prev->z.id[i];

   if(QuadInConstraint(q_prev) == false)
      return CANNOT_CONSTRAIN;
   else if( ( QuadInConstraint(q_next_copy) == false )
         && ( QuadInConstraint(q_prev) == true ) )
   {
      int RET_FLAG = ShortenVector(q_prev, q_next_copy, cVectData);

      if(RET_FLAG != CONSTR_SUCCESS) {
         ConstrVectDataReset(cVectData);
         return RET_FLAG;
      }
      for(int i = 0; i < qlen; ++i)
         q_next_copy->z.id[i] = q_prev->z.id[i] + (cVectData->tMin - BOUND_CORRECTION) * q_diff.id[i];
      COND_TEST_2;

      if(QuadInConstraintEps(q_next_copy) == true) {
         quadrature_assign(q_next_copy, q_next);
         return CONSTR_SUCCESS;
      }
      else return CONSTR_FAIL;
   }
   // if both q_prev and q_next are inside the boundary
   else {
      ConstrVectDataReset(cVectData);
      return CONSTR_NOT_NEEDED;
   }
}


// Shortens q_next vector such that every node satisfies the inequality
// A*q_next(node[i]) <=  b, provided that q_prev satisfies the constraints.
int ShortenVector(const quadrature *q_prev, const quadrature *q_next, ConstrVectData *cVectData)
{
   assert(q_prev->num_nodes == q_next->num_nodes);
   if(!QuadInConstraint(q_prev)) return CANNOT_CONSTRAIN;
   if(QuadInConstraint(q_next))  return CONSTR_NOT_NEEDED;

   int i, count;
   int k = q_prev->num_nodes;
   int dim = q_prev->dim;
   const RMatrix A = q_prev->constr->M_FULL;
   const Vector b = q_prev->constr->b_FULL;

   ConstrNodeData *cNodeDataArr = (ConstrNodeData *)malloc(k* sizeof(ConstrNodeData));
   for(int i = 0; i < k; ++i)
      cNodeDataArr[i] = ConstrNodeDataInit();

   // allocate vectors on the stack
   StaticVectorInit(z_prev_node, dim+1);
   StaticVectorInit(z_next_node, dim+1);

   for(count = 0, i = 0; i < k; ++i)
   {
      quadrature_get_elem(q_prev, i, z_prev_node);
      quadrature_get_elem(q_next, i, z_next_node);

      if( !QuadInConstraintElem(q_next, i)) {
         cNodeDataArr[i] = ShortenNode(A, b, z_prev_node, z_next_node);
         if(cNodeDataArr[i].ACTIVE == ON) ++count;
         else {
            ConstrVectDataReset(cVectData);
            free(cNodeDataArr);
            return CONSTR_FAIL;
         }
      }
   }

   if(count > 0)
   {
      Vector tVec = Vector_init(k);
      for(int i = 0; i < k; ++i)
         tVec.id[i] = cNodeDataArr[i].tMin;

      VMin vMin = VectorMin(tVec);
      int min_index             = vMin.min_index;
      cVectData->boundaryNodeId = min_index;
      cVectData->eqnId          = cNodeDataArr[min_index].eqnId;
      cVectData->tMin           = vMin.min_value;
      cVectData->ACTIVE         = ON;
      cVectData->N_OR_W         = cNodeDataArr[min_index].N_OR_W;

      Vector_free(tVec);
      free(cNodeDataArr);
      return CONSTR_SUCCESS;
   }
   else
   {
      ConstrVectDataReset(cVectData);
      free(cNodeDataArr);
      return CONSTR_UNEXPECTED;
   }
}


// Returns data information needed to map z_new onto the boundary, such that A*z_new = b_bound,
// provided that z_old satisfies A*z_old <= b_bound, and A*z_new > b_bound.
// If both A*z_new <= b_bound, and A*z_old <= b_bound,
// the routine does no further computations and default is returned.
ConstrNodeData ShortenNode(const RMatrix A, const Vector b_bound, const Vector z_old, const Vector z_new)
{
   assert(A.rows == b_bound.len);

   int i;
   int out_count;
   int nrows = A.rows;
   int ncols = A.cols;

   // Allocate vectors on stack
   StaticVectorInit(b_old, nrows);
   StaticVectorInit(b_new, nrows);
   StaticVectorInit(b_diff, nrows);
   StaticVectorInit(t, nrows);
   StaticVectorInit(dz, ncols);

   int outEqn[nrows];
   memset(outEqn, -1, nrows*sizeof(int));

   VectorAddScale(1.0, z_new, -1.0, z_old, dz);
   RMatVec(A, dz, b_diff);
   RMatVec(A, z_new, b_new);
   RMatVec(A, z_old, b_old);

   // return if z_old does not satisfy inequality constraints
   for (i = 0; i < nrows; ++i)
      if(b_old.id[i] > b_bound.id[i])
         return ConstrNodeDataInit();

   // compute t for those equations that z_new does not satisfy inequality constraints
   for (out_count = 0, i = 0; i < nrows; ++i)
      if (b_new.id[i] > b_bound.id[i]) {
         t.id[out_count] = (b_bound.id[i] - b_old.id[i]) / b_diff.id[i];
         outEqn[out_count++] = i;
      }

   // return if z_new(and z_old) satisfies constraints
   if(out_count == 0)
      return ConstrNodeDataInit();

   // most important case:
   // compute and return data if z_new does not(and z_old does) satisfy constraints
   else
   {
      int eqnId = outEqn[0];
      double tMin = t.id[0];
      for(i = 1; i < out_count; ++i)
         if(t.id[i] < tMin) {
            tMin  = t.id[i];
            eqnId = outEqn[i];
         }
      ConstrNodeData cNodeData = ConstrNodeDataInit();
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
// Most parameters are allocated on the stack due to small matrix sizes.
int ProjectNode(const CMatrix eqn_matrix, const Vector dx, Vector x_projected)
{
   int i, j, l, k;
   int M = eqn_matrix.rows;
   int N = eqn_matrix.cols;
   int LDA = M;
   int INFO;
   double TAU[MIN(M, N)];
   double WORK[N];

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

   memset(x_projected.id, 0, x_projected.len*sizeof(double));
   RMatVec(PROJECTOR, dx, x_projected);

   RMatrix_free(PROJECTOR);

   return CONSTR_SUCCESS;
}


ConstrOptData* ConstrainedOptimizationInit(quadrature *q_next)
{
   ConstrOptData *data = (ConstrOptData*)malloc(sizeof(ConstrOptData));
   data->q_next_copy = quadrature_make_full_copy(q_next);
   data->q_diff = Vector_init(q_next->z.len);
   return data;
}

void ConstrainedOptimizationFree(ConstrOptData *data)
{
   Vector_free(data->q_diff);
   quadrature_free(data->q_next_copy);
   free(data);
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

void ConstrVectDataReset(ConstrVectData *cVectData)
{
   cVectData->boundaryNodeId = -1;
   cVectData->eqnId          = -1;
   cVectData->tMin           =  1.0;
   cVectData->ACTIVE         = OFF;
   cVectData->N_OR_W         = NONE;
}

ConstrNodeData ConstrNodeDataInit()
{
   ConstrNodeData cNodeData;
   cNodeData.eqnId  = -1;
   cNodeData.tMin   =  1.0;
   cNodeData.ACTIVE = OFF;
   cNodeData.N_OR_W = NONE;

   return cNodeData;
}

__attribute__unused void ConstrNodeDataReset(ConstrNodeData *cNodeData)
{
   cNodeData->eqnId  = -1;
   cNodeData->tMin   =  1.0;
   cNodeData->ACTIVE = OFF;
   cNodeData->N_OR_W = NONE;
}


