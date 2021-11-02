#ifndef CONDITIONAL_DEBUG_H
#define CONDITIONAL_DEBUG_H


#ifdef __cplusplus
extern "C" {
#endif


#ifdef QUAD_DEBUG_ON

#define COND_TEST_1                                                                   \
double tol = POW_DOUBLE(10,-12);                                                      \
for(int i = 0; i < k; ++i)                                                            \
{                                                                                     \
   int node_index = dim * i;                                                          \
   for(int j = 0; j < num_eqns; ++j)                                                  \
   {                                                                                  \
      double test_lhs = 0.0, test_lhs_prev = 0.0;                                     \
      for(int d = 0; d < dim; ++d)                                                    \
      {                                                                               \
         test_lhs += q_next->constr->M.rid[j][d] * q_next->x[node_index+d];           \
         test_lhs_prev += q_prev->constr->M.rid[j][d] * q_prev->x[node_index+d];      \
      }                                                                               \
                                                                                      \
      if( ( fabs( test_lhs_prev - q_next->constr->b.id[j] ) <= tol )                  \
             && ( test_lhs  >= q_next->constr->b.id[j] + tol ) )                      \
      {                                                                               \
         PRINT_ERR("DID NOT PROJECT PROPERLY", __LINE__, __FILE__);                   \
         PrintInt(i, "ith node");                                                     \
         PrintNodesAndWeights(q_prev, "q_prev");                                      \
         PrintNodesAndWeights((const_quadrature *)q_next, "q_next");                  \
         RET_FLAG = CONSTR_FAILURE;                                                   \
      }                                                                               \
   }                                                                                  \
}

#define COND_TEST_2                                                       \
if(cVectData->ACTIVE == OFF)                                              \
{                                                                         \
   PRINT_ERR("CONSTRAINT FAILURE OCCURRED", __LINE__, __FILE__);          \
   PrintNodesAndWeights((const_quadrature *)q_next, "q_next");            \
   PrintNodesAndWeights(q_prev, "q_prev");                                \
}

#define COND_TEST_3                                                        \
if(QuadInConstraint( (const_quadrature *)q_next_copy ) == 0)               \
{                                                                          \
   PrintInt(cVectData->boundaryNodeId, "node_id");                         \
   PRINT_ERR("DID NOT SHORTEN SUCCESSFULLY", __LINE__, __FILE__);          \
   PrintNodesAndWeights(q_prev, "q_prev");                                 \
   PrintNodesAndWeights((const_quadrature *)q_next, "q_next");             \
}

#define COND_TEST_4                                                                   \
bool TEST_QR = TestQR(QW.rows, QW.cols, QW.id);                                       \
                                                                                      \
if(TEST_QR == FAILED)                                                                 \
{                                                                                     \
   PRINT_ERR("FAILED QR TEST, returning from NodeElimination", __LINE__, __FILE__);   \
   free(REFL);                                                                        \
   free(WORK);                                                                        \
   free(workQR);                                                                      \
   CMatrix_free(J);                                                                   \
   CMatrix_free(J_TR);                                                                \
   CMatrix_free(QW);                                                                  \
   break;                                                                             \
}

#else

#define COND_TEST_1
#define COND_TEST_2
#define COND_TEST_3
#define COND_TEST_4

#endif


#ifdef __cplusplus
}
#endif

#endif
