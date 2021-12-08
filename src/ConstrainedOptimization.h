
#ifndef CONSTRAINED_OPTIMIZATION_H
#define CONSTRAINED_OPTIMIZATION_H

#include "GENERAL_QUADRATURE.h"

#ifdef __cplusplus
extern "C" {
#endif

#define CONSTR_NOT_NEEDED 1
#define CANNOT_CONSTRAIN 2
#define CONSTR_SUCCESS  0

#define CONSTR_FAIL -1
#define CONSTR_UNEXPECTED -2

ConstrVectData ConstrVectDataInit();
void ConstrVectDataReset(ConstrVectData *cVectData);
ConstrNodeData ConstrNodeDataInit();
ATTR_UNUSED void ConstrNodeDataReset(ConstrNodeData *cNodeData);
int ConstrainedProjection(const quadrature *q_prev, quadrature *q_next);

typedef struct
{
   quadrature *q_next_copy;
   Vector q_diff;
} ConstrOptData;
ConstrOptData *ConstrainedOptimizationInit(quadrature *q_next);
void ConstrainedOptimizationFree(ConstrOptData *data);
int ConstrainedOptimization(ConstrOptData *data, const quadrature *q_prev, quadrature *q_next, ConstrVectData *cVectData);

int ShortenVector(const quadrature *q_prev, const quadrature *q_next, ConstrVectData *cVectData);
ConstrNodeData ShortenNode(const RMatrix A, const Vector b_bound, const Vector z_old, const Vector z_new);
int ProjectNode(const CMatrix eqn_matrix, const Vector dx, Vector x_projected);

#ifdef __cplusplus
}
#endif

#endif
