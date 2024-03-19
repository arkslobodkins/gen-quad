// Arkadijs Slobodkins
// Johannes Tausch
// SMU Mathematics
// Copyright 2022

// Functions that compute tensor products and helper routines
// for obtaining recursive initial guesses.

#pragma once

#if __cplusplus < 201402L
#error requires c++14 or higher
#else

#include "quadrature.hpp"

namespace gquad {

QuadCube GaussTensorCube(gq_int deg, gq_int dim);
QuadSimplex GaussTensorSimplex(gq_int deg, gq_int dim);
QuadCubeSimplex GaussTensorCubeSimplex(gq_int deg, gq_int dim1, gq_int dim2);
QuadSimplexSimplex GaussTensorSimplexSimplex(gq_int deg, gq_int dim1, gq_int dim2);

void NodesTensor2D(const QuadInterval& q1, const QuadInterval& q2, QuadDomain& q3);
void WeightsTensor2D(const QuadInterval& q1, const QuadInterval& q2, QuadDomain& q3);
void Tensor2D(const QuadInterval& q1, const QuadInterval& q2, QuadDomain& q3);
void MixedTensor(const QuadDomain& q1, const QuadDomain& q2, QuadDomain& q3);

void AddLineCube(const QuadInterval& q1, const QuadDomain& q2, QuadDomain& q3);
void AddLineSimplex(const QuadInterval& q1, const QuadDomain& q2, QuadDomain& q3);
void AddLinePyramid(const QuadInterval& q1, const QuadDomain& q2, QuadDomain& q3);

}  // namespace gquad

#endif

