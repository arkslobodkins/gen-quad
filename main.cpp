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
   // Another option is to call ReadOmega()(also in mesh.cpp) routine and supply 3 files that contain appropriate
   // arguments to constructor of Omega2D.
   if(ComputeAndOutputOmega2D(CreateIrregTrapezoid(), degrees) != 0) {
      assert(false);
   }

   if(ComputeAndOutputPentagons(degrees, SearchWidth{1}) != 0) {
      assert(false);
   }
   if(ComputeAndOutputHexagons(degrees, SearchWidth{1}) != 0) {
      assert(false);
   }

   if(ComputeAndOutputPyramids3D(degrees) != 0) {
      assert(false);
   }

   gq_int dim = 3;
   if(ComputeAndOutputCubes(degrees, dim) != 0) {
      assert(false);
   }
   if(ComputeAndOutputSimplexes(degrees, dim) != 0) {
      assert(false);
   }

   // C2 x T2
   gq_int dim1 = 2;
   gq_int dim2 = 2;
   if(ComputeAndOutputCubeSimplexes(degrees, dim1, dim2) != 0) {
      assert(false);
   }

   // T2 x T2
   if(ComputeAndOutputSimplexSimplexes(degrees, dim1, dim2) != 0) {
      assert(false);
   }

   return EXIT_SUCCESS;
}
