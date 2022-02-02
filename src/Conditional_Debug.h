#ifndef CONDITIONAL_DEBUG_H
#define CONDITIONAL_DEBUG_H

#include "Print.h"

#ifdef __cplusplus
extern "C" {
#endif


#ifdef QUAD_DEBUG_ON

#define COND_TEST_1                                                       \
for(int i = 0; i < k; ++i)                                                \
{                                                                         \
   for(int j = 0; j < num_eqns; ++j)                                      \
   if(QuadEqnOnTheBoundary(q_prev, i, j) == true                          \
       && QuadInDomainEqnElemEps(q_next, i, j) == false)                  \
      {                                                                   \
         PRINT_ERR("DID NOT PROJECT PROPERLY", __LINE__, __FILE__);       \
         PrintInt(i, "ith node");                                         \
         PrintNodeAndWeight(i, q_prev, "q_prev");                         \
         PrintNodeAndWeight(i, q_next, "q_next");                         \
      }                                                                   \
}


#define COND_TEST_2                                                           \
if(cVectData->ACTIVE == ON) {                                                 \
   if(QuadInConstraintEps(q_next_copy) == false)                              \
   {                                                                          \
      PrintInt(cVectData->boundaryNodeId, "node_id");                         \
      PRINT_ERR("DID NOT SHORTEN SUCCESSFULLY", __LINE__, __FILE__);          \
      PrintNodeAndWeight(cVectData->boundaryNodeId, q_prev, "q_prev");        \
      PrintNodeAndWeight(cVectData->boundaryNodeId, q_next_copy, "q_next");   \
   }                                                                          \
}

#define COND_TEST_4                                                            \
bool TEST_QR = TestQR(QWEIGHT);                                                \
if(TEST_QR == FAILED)                                                          \
   PRINT_ERR("FAILED QR TEST", __LINE__, __FILE__);

#define COND_PRINT_ERR(flag) \
if(flag != GQ_SUCCESS) PRINT_ERR(ERR_STRING(flag), __LINE__, __FILE__);

#define COND_TEST_5                                                           \
double __norm = 0.0;                                                          \
for(int i = 0; i < n; ++i) __norm += SQUARE(q->w[i]-q_guess->w[i]);           \
__norm = sqrt(__norm);                                                        \
if( __norm - fabs(min_value) > POW_DOUBLE(10.0, -13) )                        \
   PRINT_ERR("TEST 5: failed norm comparison", __LINE__, __FILE__);           \
if( !QuadPosWeightsEps(q_guess) )                                             \
   PRINT_ERR("TEST 5: failed positivity of weights", __LINE__, __FILE__);     \
if( fabs(V_TwoNorm(MinV) - 1.0) > POW_DOUBLE(10.0, -13) )                     \
   PRINT_ERR("TEST 5: failed orthogonality test for eigenvector", __LINE__, __FILE__);

#else // QUAD_DEBUG_OFF

#define COND_TEST_1
#define COND_TEST_2
#define COND_TEST_3
#define COND_TEST_4
#define COND_TEST_5
#define COND_PRINT_ERR

#endif


#ifdef __cplusplus
}
#endif

#endif
