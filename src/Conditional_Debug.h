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
      double m_values_test = 0.0;                                 \
         for(int d = 0; d < dim; ++d)                             \
      {                                                           \
                                                                  \
         test_lhs += cons->M.id [j] [d] * q_next->x[node_index+d];                                                 \
         test_lhs_prev += cons->M.id[j][d] * q_prev->x[node_index+d];                                              \
         m_values_test += fabs(cons->M.id[j][d]);                                                                  \
      }                                                                                                           \
                                                                                                                  \
      if( ( fabs( test_lhs_prev - cons->b.id[j] ) <= tol )  && ( test_lhs  >= cons->b.id[j] + tol ) )               \
      {                                                                                                           \
         PRINT_ERR((char *)"DID NOT PROJECT PROPERLY", 1, __LINE__, __FILE__);                                    \
         PRINT(i, "ith node");                                                                                    \
         PRINT(test_lhs, "test_lhs");                                                                             \
         PRINT(test_lhs, "test_lhs_prev");                                                                        \
         PrintNodesAndWeights(q_prev, "q_prev");                                                                  \
         PrintNodesAndWeights((const_quadrature *)q_next, "q_next");                                              \
      }                                                                                                           \
   }                                                                                                              \
}

#define COND_TEST_2                                                             \
if(constr_node_d.flag == 0)                                                     \
{                                                                               \
   PRINT_ERR((char *)"CONSTRAINT FAILURE OCCURRED", 0, __LINE__, __FILE__);     \
   PrintNodesAndWeights((const_quadrature *)q_next, "q_next");                  \
   PrintNodesAndWeights(q_prev, "q_prev");                                      \
}

#define COND_TEST_3                                                                    \
      if(InConstraint( dom_funcs, (const_quadrature *)q_next ) == 0)  \
      {                                                                                \
         PRINT(constr_node_d.node_id, "node_id");                                      \
         PRINT_ERR((char *)"DID NOT SHORTEN SUCCESSFULLY", 3, __LINE__, __FILE__);     \
         PrintNodesAndWeights(q_prev, "q_prev");                                       \
         PrintNodesAndWeights((const_quadrature *)q_next, "q_next");                   \
      }                                                                                \

#define COND_TEST_4                                                                    \
      bool TEST_QR = TestQR(ncols, weight_cols, QWEIGHT.id);                           \
                                                                                       \
      if(TEST_QR == FAILED)                                                            \
      {                                                                                \
         PRINT_ERR((char *)"FAILED QR TEST, returning", 3, __LINE__, __FILE__);        \
         break;                                                                        \
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
