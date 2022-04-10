/* Arkadijs Slobodkins
 * SMU Mathematics
 * Copyright 2022
 */


#ifndef CONSTRAINED_OPTIMIZATION_H
#define CONSTRAINED_OPTIMIZATION_H

#include "GENERAL_QUADRATURE.h"
#include "Quadrature.h"

#ifdef __cplusplus
extern "C" {
#endif

#define CONSTR_NOT_NEEDED 1
#define CANNOT_CONSTRAIN 2
#define CONSTR_SUCCESS  0

#define CONSTR_FAIL -1
#define CONSTR_UNEXPECTED -2

typedef enum { NODE, WEIGHT, NONE } NODE_OR_WEIGHT;

typedef struct
{
   int            eqnId;
   double         tMin;
   bool_enum      ACTIVE;
   NODE_OR_WEIGHT N_OR_W;
} ConstrNodeData;

typedef struct
{
   int            boundaryNodeId;
   int            eqnId;
   double         tMin;
   bool_enum      ACTIVE;
   NODE_OR_WEIGHT N_OR_W;
} ConstrVectData;

typedef struct
{
   quadrature *q_next_copy;
   Vector q_diff;
} ConstrOptData;


ConstrVectData ConstrVectDataInit();
void ConstrVectDataReset(ConstrVectData *cVectData);
ConstrNodeData ConstrNodeDataInit();
__attribute__unused void ConstrNodeDataReset(ConstrNodeData *cNodeData);


int ConstrainedProjection(const quadrature *q_prev, quadrature *q_next);


ConstrOptData *ConstrainedOptimizationInit(quadrature *q_next);
void ConstrainedOptimizationFree(ConstrOptData *data);
int ConstrainedOptimization(ConstrOptData *data, const quadrature *q_prev, quadrature *q_next, ConstrVectData *cVectData);


// Shortens q_next vector such that every node satisfies the inequality
// A*q_next(node[i]) <=  b, provided that q_prev satisfies the constraints.
int ShortenVector(const quadrature *q_prev, const quadrature *q_next, ConstrVectData *cVectData);

// Returns data information needed to map z_new onto the boundary, such that A*z_new = b_bound,
// provided that z_old satisfies A*z_old <= b_bound, and A*z_new > b_bound.
// If both A*z_new <= b_bound, and A*z_old <= b_bound,
// the routine does no further computations and default is returned.
ConstrNodeData ShortenNode(const RMatrix A, const Vector b_bound, const Vector z_old, const Vector z_new);

// Projects dx onto equations specified by eqn_matrix. Reduced matrix Q of eqn_matrix
// is extracted from QR factorization, and projector P is computed by
// P = I-Q_reduced*Q_reduced'. Results are stored in x_projected.
// Most parameters are allocated on the stack due to small matrix sizes.
int ProjectNode(const CMatrix eqn_matrix, const Vector dx, Vector x_projected);

#ifdef __cplusplus
}
#endif

#endif
