#ifndef SETCONSTRAINT_H
#define SETCONSTRAINT_H

#include "GENERAL_QUADRATURE.h"

#ifdef __cplusplus
extern "C" {
#endif

void SetConstraintCube(ConstraintFuncs *constr);
void SetConstraintSimplex(ConstraintFuncs *constr);
void SetConstraintCubeSimplex(ConstraintFuncs *constr);
void SetConstraintSimplexSimplex(ConstraintFuncs *constr);
void SetConstraintCubeSimplexSimplex(ConstraintFuncs *constr);

#ifdef __cplusplus
}
#endif

#endif
