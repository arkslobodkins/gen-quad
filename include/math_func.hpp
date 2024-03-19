// Arkadijs Slobodkins
// Johannes Tausch
// SMU Mathematics
// Copyright 2022

#pragma once

#if __cplusplus < 201402L
#error requires c++14 or higher
#else

#include <omp.h>

#include <cmath>

namespace gquad { namespace math {

// fast power function
// may introduce slighly greater roundoff
// error than the standard pow function
template <typename RealType>
inline RealType IntegerPower(RealType x, long int power) {
   RealType res{1};
   for(;;) {
      if(power & 1) {
         res *= x;
      }
      power >>= 1;
      if(!power) {
         break;
      }
      x *= x;
   }
   return res;
}

inline long int factorial(long int n) {
   return n < 2L ? 1L : n * factorial(n - 1L);
}

inline long int binomial(long int k, long int n) {
   return factorial(n) / (factorial(k) * factorial(n - k));
}

// evaluate exponential product of doubles in long double precision
inline long double expNDim(long int dim, const double x[]) {
   long double expVal = 1.L;
   for(long int i = 0; i < dim; ++i) {
      expVal *= std::exp(
          (long double)x[i]);  // expl is not used since it might not be supported on older compilers
   }
   return expVal;
}

inline long double expIntegral1D(long double c) {
   return 1.L / c * (std::exp(c) - 1.L);
}

}}  // namespace gquad::math

#endif

