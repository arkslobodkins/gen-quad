// Arkadijs Slobodkins
// Johannes Tausch
// SMU Mathematics
// Copyright 2022

#include <cassert>
#include <cstdlib>
#include <mesh.hpp>
#include <quad_driver.hpp>

using namespace gquad;

int main() {
   StdVector<gq_int> degrees = sequence(5, 7);  // degrees 5, 6, 7

   // examples:
   // The first argument to ComputeAndOutputOmega2D is an object of type Omega2D,
   // which is defined in include/mesh.hpp. CreateIrregTrapezoid(), CreateIrreg5(), CreatePentagon(),
   // and CreateHexagon() defined in src/mesh.cpp demonstrate how to construct such objects.
   // Another option is to call ReadOmega()(also in mesh.cpp) routine and supply 3 files that contain
   // appropriate arguments to constructor of Omega2D.

   assert(ComputeAndOutputIntervals(degrees) == 0);
//   assert(ComputeAndOutputOmega2D(CreateIrreg5(), degrees, SearchWidth{1}) == 0);
//   assert(ComputeAndOutputPentagons(degrees, SearchWidth{1}) == 0);
//   assert(ComputeAndOutputHexagons(degrees, SearchWidth{1}) == 0);
//
//   gq_int dim = 3;
//   assert(ComputeAndOutputPyramids3D(degrees) == 0);
//   assert(ComputeAndOutputCubes(degrees, dim) == 0);
//   assert(ComputeAndOutputSimplexes(degrees, dim) == 0);
//
//   gq_int dim1 = 2;
//   gq_int dim2 = 2;
//   assert(ComputeAndOutputCubeSimplexes(degrees, dim1, dim2) == 0);     // C2 x T2
//   assert(ComputeAndOutputSimplexSimplexes(degrees, dim1, dim2) == 0);  // T2 x T2
//
   return EXIT_SUCCESS;
}
