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
   StdVector<gq_int> degrees = sequence(43, 50, 2);  // degrees 5, 6, 7

   // examples:
   // The first argument to ComputeAndOutputOmega2D is an object of type Omega2D,
   // which is defined in include/mesh.hpp. CreateIrregTrapezoid(), CreateIrreg5(), CreatePentagon(),
   // and CreateHexagon() defined in src/mesh.cpp demonstrate how to construct such objects.
   // Another option is to call ReadOmega()(also in mesh.cpp) routine and supply 3 files that contain appropriate
   // arguments to constructor of Omega2D.
   if(ComputeAndOutputPentagons(degrees) != 0) {
      assert(false);
   }

   return EXIT_SUCCESS;
}
