/* Arkadijs Slobodkins
 * SMU Mathematics
 * Copyright 2022
 */


#ifndef LEAST_SQUARES_NEWTONPLASMA_H
#define LEAST_SQUARES_NEWTONPLASMA_H

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


// LeastSquaresNewton
// Receives initial quadrature guess. Primarily solves
// underdetermined systems of equations in the least
// squares sense. Returns success if algorithm converged
// and all nodes are inside of the domain and if all nodes are positive.
LSQ_out LeastSquaresNewton(const bool_enum CONSTR_OPT, quadrature *q_orig);
LSQ_out LevenbergMarquardt(const bool_enum CONSTR_OPT, quadrature *q_orig);

#ifdef __cplusplus
}
#endif

#endif


