// Arkadijs Slobodkins
// Johannes Tausch
// SMU Mathematics
// Copyright 2022

// Routines that generate recursive initial guesses and run Node Elimination algorithm.
// Quadrature rules and Node Elimination History for all calls to Node Elimination are returned.

#pragma once

#if __cplusplus < 201402L
#error requires c++14 or higher
#else

#include "quadrature.hpp"

namespace gquad {


// Nothing to be done, return Gauss-Legendre quadrature on [0, 1].
QuadInterval ComputeInterval(gq_int deg);


// computes 3-D pyramid
std::pair<QuadPyramid3D, StdVector<History>> ComputePyramid3D(gq_int deg, gq_int search_width);


// Routines below run Node Elimination algorithm using recursive initial guess.
std::pair<QuadCube, StdVector<History>> ComputeCube(gq_int deg, gq_int dim, gq_int search_width);
std::pair<QuadSimplex, StdVector<History>> ComputeSimplex(gq_int deg, gq_int dim, gq_int search_width);
std::pair<QuadCubeSimplex, StdVector<History>> ComputeCubeSimplex(gq_int deg, gq_int dim1, gq_int dim2,
                                                                  gq_int search_width);


// Runs Node Elimination algorithm using recursive initial guess for
// dim-1 simplex and dim-2 simplex, but takes a tensor product of simplex1 x simplex2.
std::pair<QuadSimplexSimplex, StdVector<History>> ComputeSimplexSimplex(gq_int deg, gq_int dim1, gq_int dim2,
                                                                        gq_int search_width);


std::pair<QuadOmega2D, StdVector<History>> ComputeOmega2D(Omega2D omega, gq_int deg, gq_int search_width);
std::pair<QuadOmega2D, StdVector<History>> ComputePentagon(gq_int deg, gq_int search_width);
std::pair<QuadOmega2D, StdVector<History>> ComputeHexagon(gq_int deg, gq_int search_width);


}  // namespace gquad

#endif

