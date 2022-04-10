/* Arkadijs Slobodkins
 * SMU Mathematics
 * Copyright 2022
 */


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
   int k = q_prev->num_nodes;
   int dim = q_prev->dim;
   assert(dim == q_next->constr->M.cols);

   int num_eqns = q_next->constr->M.rows;
   RMatrix M = q_next->constr->M;

   Tensor2D x_next_T = DoubleToTensor2D(k, dim, q_next->x);
   Tensor2D x_prev_T = DoubleToTensor2D(k, dim, q_prev->x);

   StaticVectorInit(dim, node_change);
   StaticVectorInit(dim, node_projected);
   Tensor1D node_change_T = VectorToTensor1D(node_change);
   Tensor1D node_projected_T = VectorToTensor1D(node_projected);

   for(int i = 0; i < k; ++i)
   {
      int active_eqn_count = 0;
      bool do_project = false;
      bool eqn_flags[num_eqns];

      for(int j = 0; j < num_eqns; ++j) eqn_flags[j] = false;

      if( QuadOnTheBoundary(q_prev, i)  && !QuadInDomainElem(q_next, i) )
         do_project = true;

      if(do_project)
      {
         for(int j = 0; j < num_eqns; ++j)
         {
            if(QuadEqnOnTheBoundary(q_prev, i, j))
            {
               ++active_eqn_count;
               eqn_flags[j] = true;
            }
         }
         CMatrix eqn_matrix = CMatrix_init(dim, active_eqn_count);
         int count = 0;
         for(int j = 0; j < num_eqns; ++j)
         {
            if(eqn_flags[j] == true)
            {
               for(int d = 0; d < dim; ++d)
                  eqn_matrix.cid[count][d]= M.rid[j][d];
               ++count;
            }
         }

         for(int d = 0; d < dim; ++d)
            TID1(node_change_T, 1) = TID2(x_next_T, i, d) - TID2(x_prev_T, i, d);
         int P_FLAG = ProjectNode(eqn_matrix, node_change, node_projected);
         if(P_FLAG != CONSTR_SUCCESS)
         {
            CMatrix_free(eqn_matrix);
            return P_FLAG;
         }
         for(int d = 0; d < dim; ++d)
            TID2(x_next_T, i, d) = TID2(x_prev_T, i, d) + TID1(node_projected_T, d) - POW_DOUBLE(10, -13);

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

   Vector q_diff = data->q_diff;
   quadrature *q_next_copy = data->q_next_copy;
   Vector x_next_copy = q_next_copy->z;
   Vector x_prev = q_prev->z;

   for(int i = 0; i < qlen; ++i)
      q_diff.id[i] = x_next_copy.id[i] - x_prev.id[i];

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
         x_next_copy.id[i] = x_prev.id[i] + (cVectData->tMin - BOUND_CORRECTION) * q_diff.id[i];
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


int ShortenVector(const quadrature *q_prev, const quadrature *q_next, ConstrVectData *cVectData)
{
   assert(q_prev->num_nodes == q_next->num_nodes);
   if(!QuadInConstraint(q_prev)) return CANNOT_CONSTRAIN;
   if(QuadInConstraint(q_next))  return CONSTR_NOT_NEEDED;

   int k = q_prev->num_nodes;
   int dim = q_prev->dim;
   const RMatrix A = q_prev->constr->M_FULL;
   const Vector b = q_prev->constr->b_FULL;

   ConstrNodeData *cNodeDataArr = (ConstrNodeData *)malloc(k* sizeof(ConstrNodeData));
   for(int i = 0; i < k; ++i)
      cNodeDataArr[i] = ConstrNodeDataInit();

   // allocate vectors on the stack
   StaticVectorInit(dim+1, z_prev_node);
   StaticVectorInit(dim+1, z_next_node);

   int count = 0;
   for(int i = 0; i < k; ++i)
   {
      quadrature_get_elem(q_prev, i, z_prev_node);
      quadrature_get_elem(q_next, i, z_next_node);

      if( !QuadInConstraintElem(q_next, i) ) {
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


ConstrNodeData ShortenNode(const RMatrix A, const Vector b_bound, const Vector z_old, const Vector z_new)
{
   assert(A.rows == b_bound.len);

   int i;
   int out_count;
   int nrows = A.rows;
   int ncols = A.cols;

   // Allocate vectors on stack
   StaticVectorInit(nrows, b_old);
   StaticVectorInit(nrows, b_new);
   StaticVectorInit(nrows, b_diff);
   StaticVectorInit(nrows, t);
   StaticVectorInit(ncols, dz);

   int outEqn[nrows];
   memset(outEqn, -1, nrows*sizeof(int));

   VAddScale(1.0, z_new, -1.0, z_old, dz);
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


int ProjectNode(const CMatrix eqn_matrix, const Vector dx, Vector x_projected)
{
   int M = eqn_matrix.rows;
   int N = eqn_matrix.cols;

   StaticVectorInit(MIN(M, N), TAU);

   int INFO = DGEQR2_LAPACK(eqn_matrix, TAU);
   if(INFO != 0) return CONSTR_FAIL;

   INFO = DORGQR_LAPACK(eqn_matrix, TAU);
   if(INFO != 0) return CONSTR_FAIL;

   /******************************************
   \* Obtain P = I - Q_REDUCED x Q_REDUCED^ \*
   ******************************************/

   double Q_REDUCED[M][N];
   for(int i = 0; i < M; ++i)
      for(int j = 0; j < N; ++j)
         Q_REDUCED[i][j] = eqn_matrix.cid[j][i];

   double Q_TRANS[N][M];
   for(int i = 0; i < M; ++i)
      for(int j = 0; j < N; ++j)
         Q_TRANS[j][i] = Q_REDUCED[i][j];

   RMatrix PROJECTOR = RMatrix_init(M, M);
   for(int i = 0; i < M; ++i)
      for(int k = 0; k < N; ++k)
         for(int j = 0; j < M; ++j)
            PROJECTOR.rid[i][j] += Q_REDUCED[i][k]*Q_TRANS[k][j];

   for(int i = 0; i < M; ++i)
      for(int j = 0; j < M; ++j)
         PROJECTOR.rid[i][j] = -PROJECTOR.rid[i][j];
   for(int i = 0; i < M; ++i) ++PROJECTOR.rid[i][i];

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
   cVectData.tMin           = 1.0;
   cVectData.ACTIVE         = OFF;
   cVectData.N_OR_W         = NONE;

   return cVectData;
}

void ConstrVectDataReset(ConstrVectData *cVectData)
{
   cVectData->boundaryNodeId = -1;
   cVectData->eqnId          = -1;
   cVectData->tMin           = 1.0;
   cVectData->ACTIVE         = OFF;
   cVectData->N_OR_W         = NONE;
}

ConstrNodeData ConstrNodeDataInit()
{
   ConstrNodeData cNodeData;
   cNodeData.eqnId  = -1;
   cNodeData.tMin   = 1.0;
   cNodeData.ACTIVE = OFF;
   cNodeData.N_OR_W = NONE;

   return cNodeData;
}

__attribute__unused void ConstrNodeDataReset(ConstrNodeData *cNodeData)
{
   cNodeData->eqnId  = -1;
   cNodeData->tMin   = 1.0;
   cNodeData->ACTIVE = OFF;
   cNodeData->N_OR_W = NONE;
}


