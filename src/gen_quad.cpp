// Arkadijs Slobodkins
// Johannes Tausch
// SMU Mathematics
// Copyright 2022

#if __cplusplus < 201402L
#error requires c++14 or higher
#else

#include "../include/gen_quad.hpp"

namespace gquad {

StdVector<gq_int> sequence(gq_int first, gq_int last, gq_int stride) {
   GEN_QUAD_ASSERT_ALWAYS(last >= first);
   GEN_QUAD_ASSERT_ALWAYS(stride >= 1);
   StdVector<gq_int> seq{};
   for(gq_int i = first; i <= last; i += stride) {
      seq.push_back(i);
   }
   return seq;
}

bool is_reasonable_degree_and_dimension(gq_int deg, gq_int dim) {
   GEN_QUAD_ASSERT_ALWAYS(deg > 0);
   GEN_QUAD_ASSERT_ALWAYS(dim > 0);

   if(dim == 1 && deg <= 100) {
      return true;
   }
   if(dim == 2 && deg <= 100) {
      return true;
   } else if(dim == 3 && deg <= 50) {
      return true;
   } else if(dim == 4 && deg <= 25) {
      return true;
   } else if(dim == 5 && deg <= 20) {
      return true;
   } else if(dim == 6 && deg <= 15) {
      return true;
   } else if(dim == 7 && deg <= 10) {
      return true;
   } else if(dim == 8 && deg <= 8) {
      return true;
   } else {
      return false;
   }
}

#ifdef _OPENMP
bool omp_condition(gq_int deg, gq_int dim) {
   if(dim == 2 && deg >= 17) {
      return true;
   } else if(dim == 3 && deg >= 11) {
      return true;
   } else if(dim == 4 && deg >= 7) {
      return true;
   } else if(dim == 5 && deg >= 6) {
      return true;
   } else if(dim == 6 && deg >= 5) {
      return true;
   } else if(dim > 6) {
      return true;
   } else {
      return false;
   }
}
#endif

}  // namespace gquad

#endif

