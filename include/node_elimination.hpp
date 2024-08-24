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


// Node Elimination driver routine. It expects input to be initial(preferably exact)
// quadrature rule of the desired degree and dimension.
std::pair<std::unique_ptr<QuadDomain>, History> NodeElimination(const QuadDomain& quad_init,
                                                                gq_int search_width);


struct PredictorTimer {
   static double predictor_time_total;
};


}  // namespace gquad

#endif

