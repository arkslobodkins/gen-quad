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

struct Penalty {
   double factor;
   double distance;
};

LsqOut LeastSquaresNewton(QuadDomain& q, Basis& basis);
Penalty PenaltyGradientFac(gq_int skip, const QuadIdealPolytope& q, double dzFac, const QuadArray& dz,
                           const QuadArray& dg);
double PenaltyGradientFac(const QuadOmega2D& q, const QuadArray& dz, const QuadArray& dg);

QuadArray ComputePenaltyGradient(const QuadIdealPolytope& q);
QuadArray ComputePenaltyGradient(const QuadOmega2D& q, bool do_inner);

struct LsqSolveTimer {
   static double lsq_time_total;
};

}  // namespace gquad

#endif

