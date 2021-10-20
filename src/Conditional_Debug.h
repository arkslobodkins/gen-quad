#ifndef CONDITIONAL_DEBUG_H
#define CONDITIONAL_DEBUG_H


#ifdef __cplusplus
extern "C" {
#endif


#ifdef QUAD_DEBUG_ON

#define COND_TEST_1                                               \
for(int i = 0; i < k; ++i)                                        \
{                                                                 \
   int node_index = dim * i;                                      \
   for(int j = 0; j < num_eqns; ++j)                              \
   {                                                              \
      double test_lhs = 0.0, test_lhs_prev = 0.0;                 \
      for(int d = 0; d < dim; ++d)                                \
      {                                                           \
         test_lhs += q_next->cons->M.id [j] [d] * q_next->x[node_index+d];                                             \
         test_lhs_prev += q_prev->cons->M.id[j][d] * q_prev->x[node_index+d];                                          \
      }                                                                                                                \
                                                                                                                       \
      if( ( fabs( test_lhs_prev - q_next->cons->b.id[j] ) <= tol )  && ( test_lhs  >= q_next->cons->b.id[j] + tol ) )  \
      {                                                                                                                \
         PRINT_ERR("DID NOT PROJECT PROPERLY", __LINE__, __FILE__);                                                    \
         PrintInt(i, "ith node");                                                                                      \
         PrintNodesAndWeights(q_prev, "q_prev");                                                                       \
         PrintNodesAndWeights((const_quadrature *)q_next, "q_next");                                                   \
      }                                                                                                                \
   }                                                                                                                   \
}

#define COND_TEST_2                                                              \
if(cNodeData.flag == 0)                                                          \
{                                                                                \
   PRINT_ERR("CONSTRAINT FAILURE OCCURRED", __LINE__, __FILE__);                 \
   PrintNodesAndWeights((const_quadrature *)q_next, "q_next");                   \
   PrintNodesAndWeights(q_prev, "q_prev");                                       \
}

#define COND_TEST_3                                                              \
if(q_next->inConstraint( (const_quadrature *)q_next ) == 0)                      \
{                                                                                \
   PrintInt(cNodeData.nodeId, "node_id");                                        \
   PRINT_ERR("DID NOT SHORTEN SUCCESSFULLY", __LINE__, __FILE__);                \
   PrintNodesAndWeights(q_prev, "q_prev");                                       \
   PrintNodesAndWeights((const_quadrature *)q_next, "q_next");                   \
}

#define COND_TEST_4                                                                   \
bool TEST_QR = TestQR(ncols, weight_cols, QWEIGHT.id);                                \
                                                                                      \
if(TEST_QR == FAILED)                                                                 \
{                                                                                     \
   PRINT_ERR("FAILED QR TEST, returning from NodeElimination", __LINE__, __FILE__);   \
   free(TAU);                                                                         \
   free(WORK);                                                                        \
   free(workQR);                                                                      \
   Vector_free(JACOBIAN);                                                             \
   Vector_free(QWEIGHT);                                                              \
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
