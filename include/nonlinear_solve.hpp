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


enum class SolFlag { found, not_found };


struct LsqOut {
   SolFlag sol_flag;
   gq_int its;
   double residual;
};


LsqOut LeastSquaresNewton(QuadDomain& q, Basis& basis);


struct LsqSolveTimer {
   static double lsq_time_total;
};


}  // namespace gquad

#endif

