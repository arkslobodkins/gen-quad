/* Arkadijs Slobodkins
 * SMU Mathematics
 * Copyright 2022
 */


#ifndef CONDITIONAL_DEBUG_H
#define CONDITIONAL_DEBUG_H

#include "Print.h"

#ifdef __cplusplus
extern "C" {
#endif


#ifdef QUAD_DEBUG_ON

#define COND_TEST_1                                                        \
for(int i = 0; i < k; ++i)                                                 \
{                                                                          \
   for(int j = 0; j < num_eqns; ++j)                                       \
   {                                                                       \
      if(QuadEqnOnTheBoundary(q_prev, i, j) == true                        \
         && QuadInDomainEqnElemEps(q_next, i, j) == false)                 \
      {                                                                    \
         PRINT_ERR("DID NOT PROJECT PROPERLY", __LINE__, __FILE__);        \
         PrintInt(i, "ith node");                                          \
         PrintNodeAndWeight(i, q_prev, "q_prev");                          \
         PrintNodeAndWeight(i, q_next, "q_next");                          \
      }                                                                    \
   }                                                                       \
}


#define COND_TEST_2                                                           \
if(cVectData->ACTIVE == ON)                                                   \
{                                                                             \
   if(QuadInConstraintEps(q_next_copy) == false)                              \
   {                                                                          \
      PrintInt(cVectData->boundaryNodeId, "node_id");                         \
      PRINT_ERR("DID NOT SHORTEN SUCCESSFULLY", __LINE__, __FILE__);          \
      PrintNodeAndWeight(cVectData->boundaryNodeId, q_prev, "q_prev");        \
      PrintNodeAndWeight(cVectData->boundaryNodeId, q_next_copy, "q_next");   \
   }                                                                          \
}


#define COND_TEST_3                                                           \
bool TEST_QR = TestQR(QFull);                                                 \
if(TEST_QR == QR_FAILED)                                                      \
   PRINT_ERR("QR_FAILED QR TEST", __LINE__, __FILE__);


#define COND_TEST_4                                                                        \
for(int i = 0; i < Z.rows; ++i)                                                            \
   if( fabs(Z.rid[i][whichIndex[i]]) > POW_DOUBLE(10, -12) )                               \
      PRINT_ERR("TEST 4: failed positivity of weights", __LINE__, __FILE__);               \
                                                                                           \
for(int i = 0; i < VT.rows; ++i)                                                           \
{                                                                                          \
   double __norm = 0.0;                                                                    \
   for(int j = 0; j < n_cur; ++j)                                                          \
      __norm += SQUARE(VT.cid[j][i]);                                                      \
   if( fabs(__norm - 1.0) > POW_DOUBLE(10.0, -13) )                                        \
      PRINT_ERR("TEST 4: failed orthogonality test for eigenvectors", __LINE__, __FILE__); \
}


#define COND_PRINT_ERR(flag) \
if(flag != GQ_SUCCESS) PRINT_ERR(ERR_STRING(flag), __LINE__, __FILE__);

#else // QUAD_DEBUG_OFF

#define COND_TEST_1
#define COND_TEST_2
#define COND_TEST_3
#define COND_TEST_4
#define COND_PRINT_ERR

#endif


#ifdef __cplusplus
}
#endif

#endif
