// Arkadijs Slobodkins
// Johannes Tausch
// SMU Mathematics
// Copyright 2022

#pragma once

#if __cplusplus < 201402L
#error requires c++14 or higher
#else

#include "quadrature.hpp"

namespace gquad {


// <ComputeAndOutput> prefix indicates higher level routines that compute quadrature
// rules for a vector of degrees rather than one degree at a time. These routines also
// output quadrature rules and other information to files located in results directory.

// gq_int indicates whether all quadrature rules have been computed successfully.
// 0 if succeeded, 1 if computation has been terminated due to illegal input, exceptions, etc.
gq_int ComputeAndOutputIntervals(const StdVector<gq_int>& deg);


gq_int ComputeAndOutputPyramids3D(const StdVector<gq_int>& deg, SearchWidth search_width = SearchWidth{1});


gq_int ComputeAndOutputCubes(const StdVector<gq_int>& deg, gq_int dim,
                             SearchWidth search_width = SearchWidth{1});


gq_int ComputeAndOutputSimplexes(const StdVector<gq_int>& deg, gq_int dim,
                                 SearchWidth search_width = SearchWidth{1});


gq_int ComputeAndOutputCubeSimplexes(const StdVector<gq_int>& deg, gq_int dim1, gq_int dim2,
                                     SearchWidth search_width = SearchWidth{1});


gq_int ComputeAndOutputSimplexSimplexes(const StdVector<gq_int>& deg, gq_int dim1, gq_int dim2,
                                        SearchWidth search_width = SearchWidth{1});


gq_int ComputeAndOutputPentagons(const StdVector<gq_int>& deg, SearchWidth search_width = SearchWidth{1});


gq_int ComputeAndOutputHexagons(const StdVector<gq_int>& deg, SearchWidth search_width = SearchWidth{1});


gq_int ComputeAndOutputOmega2D(Omega2D omega, const StdVector<gq_int>& deg,
                               SearchWidth search_width = SearchWidth{1});


// Returns null pointer if computation has been terminated due to illegal input, exceptions, etc.
std::unique_ptr<QuadInterval> QuadratureInterval(gq_int deg);


std::unique_ptr<QuadPyramid3D> QuadraturePyramid3D(gq_int deg, SearchWidth search_width = SearchWidth{1});


std::unique_ptr<QuadCube> QuadratureCube(gq_int deg, gq_int dim, SearchWidth search_width = SearchWidth{1});


std::unique_ptr<QuadSimplex> QuadratureSimplex(gq_int deg, gq_int dim,
                                               SearchWidth search_width = SearchWidth{1});


std::unique_ptr<QuadCubeSimplex> QuadratureCubeSimplex(gq_int deg, gq_int dim1, gq_int dim2,
                                                       SearchWidth search_width = SearchWidth{1});


std::unique_ptr<QuadSimplexSimplex> QuadratureSimplexSimplex(gq_int deg, gq_int dim1, gq_int dim2,
                                                             SearchWidth search_width = SearchWidth{1});


std::unique_ptr<QuadOmega2D> QuadraturePentagon(gq_int deg, SearchWidth search_width = SearchWidth{1});


std::unique_ptr<QuadOmega2D> QuadratureHexagon(gq_int deg, SearchWidth search_width = SearchWidth{1});


std::unique_ptr<QuadOmega2D> QuadratureOmega2D(Omega2D omega, gq_int deg,
                                               SearchWidth search_width = SearchWidth{1});


// Generates file in results/quad_rules directory.
// Useful for the second type of driver routines,
// i.e. the ones with <Quadrature> prefix.
void QuadToFile(const QuadDomain& q);


}  // namespace gquad

#endif

