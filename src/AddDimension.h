/* Arkadijs Slobodkins
 * SMU Mathematics
 * Copyright 2022
 */


#ifndef ADD_DIMENSION_H
#define ADD_DIMENSION_H

#include "GENERAL_QUADRATURE.h"
#include "Quadrature.h"

#ifdef __cplusplus
extern "C" {
#endif


// Computes tensor product of arbitrary quadrature over (dim-1)-dimensional
// domain Ω and 1-dimensional interval [0, 1]. Quadrature of dimension dim
// and degree p over interval Ω x [0,1] is generated. Stores interval quadrature nodes
// as the first coordinate.
void AddLineFirst(const quadrature *q1D, const quadrature *quad_prev, quadrature *quad_new);


// Computes tensor product of unit (dim-1)-dimensional simplex
// and 1-dimensional interval [0, 1], and maps the product to the unit
// simplex of dimension dim using Duffy transformation.
void AddLineSimplex(const quadrature *q1D, const quadrature *quad_prev, quadrature *quad_new);


// Maps nodes and weights from cube of dimension dim to
// simplex of dimension dim using generalized Duffy transformation.
void GeneralDuffy(quadrature *q);



// Computes tensor product of dim1-dimensional quadrature q1
// and dim2-dimensional quadrature q2.
void MixedTensor(const quadrature *q1, const quadrature *q2, quadrature *q_tp);


#ifdef __cplusplus
}
#endif

#endif
