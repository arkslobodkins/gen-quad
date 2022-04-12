/* Arkadijs Slobodkins
 * SMU Mathematics
 * Copyright 2022
 */


#ifndef NONLINEAR_SOLVE_H
#define NONLINEAR_SOLVE_H

#include "GENERAL_QUADRATURE.h"
#include "Quadrature.h"

#ifdef __cplusplus
extern "C" {
#endif

#define SOL_FOUND true
#define SOL_NOT_FOUND false

typedef struct
{
   int its;
   bool SOL_FLAG;
} LSQ_out;

// LeastSquaresNewton && LevenbergMarquardt
// Both solvers receive initial quadrature guess and parameter for
// enabling/disabling constrained optimization mode.
// Primarily designed for solving underdetermined systems of equations.
// Returns SOL_FOUND if algorithm converged and all nodes are inside
// of the domain and if all nodes are positive, otherwise SOL_NOT_FOUND.

LSQ_out LeastSquaresNewton(const bool_enum CONSTR_OPT, quadrature *q_orig);
LSQ_out LevenbergMarquardt(const bool_enum CONSTR_OPT, quadrature *q_orig);

#ifdef __cplusplus
}
#endif

#endif


