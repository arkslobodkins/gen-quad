// Arkadijs Slobodkins
// Johannes Tausch
// SMU Mathematics
// Copyright 2022

#if __cplusplus < 201402L
#error requires c++14 or higher
#else

#include "../include/compute_quad.hpp"

#include "../include/node_elimination.hpp"
#include "../include/quad_tensor.hpp"

namespace gquad {


QuadInterval ComputeInterval(gq_int deg) {
   GEN_QUAD_ASSERT_ALWAYS(deg > 0);
   return quadrature_gauss_legendre(deg);
}


std::pair<QuadCube, StdVector<History>> ComputeCube(gq_int deg, gq_int dim, gq_int search_width) {
   GEN_QUAD_ASSERT_ALWAYS(deg > 0);
   GEN_QUAD_ASSERT_ALWAYS(dim > 1);

   // 2-dim cube
   QuadInterval qi = quadrature_gauss_legendre(deg);
   QuadCube qc_init(deg, 2, mult_nodes(qi, qi));
   AddLineCube(qi, qi, qc_init);
   auto cur = NodeElimination(qc_init, search_width);
   auto &qc_cur = dynamic_cast<QuadCube &>(*cur.first);
   StdVector<History> hist;
   hist.push_back(cur.second);

   // d-dim cube
   for(gq_int d = 3; d <= dim; ++d) {
      QuadCube qc_init(deg, d, mult_nodes(qi, qc_cur));
      AddLineCube(qi, qc_cur, qc_init);

      auto update = NodeElimination(qc_init, search_width);
      const auto &qc_update = dynamic_cast<QuadCube &>(*update.first);
      qc_cur.reinit_copy(qc_update);
      hist.push_back(update.second);
   }

   return {qc_cur, hist};
}


std::pair<QuadSimplex, StdVector<History>> ComputeSimplex(gq_int deg, gq_int dim, gq_int search_width) {
   GEN_QUAD_ASSERT_ALWAYS(deg > 0);
   GEN_QUAD_ASSERT_ALWAYS(dim > 1);

   // 2-dim simplex
   QuadInterval qgl = quadrature_gauss_legendre(deg);
   QuadInterval qgj = quadrature_gauss_jacobi(deg, 0, 1);
   QuadSimplex qs_init(deg, 2, qgl.num_nodes() * qgj.num_nodes());
   AddLineSimplex(qgj, qgl, qs_init);

   auto cur = NodeElimination(qs_init, search_width);
   auto &qs_cur = dynamic_cast<QuadSimplex &>(*cur.first);
   StdVector<History> hist;
   hist.push_back(cur.second);

   // d-dim simplex
   for(gq_int d = 3; d <= dim; ++d) {
      QuadInterval qgj_next = quadrature_gauss_jacobi(deg, 0, d - 1);
      QuadSimplex qs_init_next(deg, d, mult_nodes(qgj_next, qs_cur));
      AddLineSimplex(qgj_next, qs_cur, qs_init_next);

      auto update = NodeElimination(qs_init_next, search_width);
      const auto &qs_update = dynamic_cast<QuadSimplex &>(*update.first);
      qs_cur.reinit_copy(qs_update);
      hist.push_back(update.second);
   }

   return {qs_cur, hist};
}


std::pair<QuadCubeSimplex, StdVector<History>> ComputeCubeSimplex(gq_int deg, gq_int dim1, gq_int dim2,
                                                                  gq_int search_width) {
   GEN_QUAD_ASSERT_ALWAYS(deg > 0);
   GEN_QUAD_ASSERT_ALWAYS(dim1 > 0);
   GEN_QUAD_ASSERT_ALWAYS(dim2 > 1);

   auto qs_h1 = ComputeSimplex(deg, dim2, search_width);
   const auto &qs = qs_h1.first;
   auto &hist = qs_h1.second;

   QuadInterval qi = quadrature_gauss_legendre(deg);
   QuadCubeSimplex qcs_init(deg, 1, dim2, mult_nodes(qi, qs));
   AddLineCube(qi, qs, qcs_init);
   auto cur = NodeElimination(qcs_init, search_width);
   auto &qcs_cur = dynamic_cast<QuadCubeSimplex &>(*cur.first);
   hist.push_back(cur.second);

   for(gq_int d = 2; d <= dim1; ++d) {
      QuadCubeSimplex qcs_init_next(deg, d, dim2, mult_nodes(qi, qcs_cur));
      AddLineCube(qi, qcs_cur, qcs_init_next);

      auto update = NodeElimination(qcs_init_next, search_width);
      const auto &qcs_update = dynamic_cast<QuadCubeSimplex &>(*update.first);
      qcs_cur.reinit_copy(qcs_update);
      hist.push_back(update.second);
   }

   return {qcs_cur, hist};
}


std::pair<QuadSimplexSimplex, StdVector<History>> ComputeSimplexSimplex(gq_int deg, gq_int dim1, gq_int dim2,
                                                                        gq_int search_width) {
   GEN_QUAD_ASSERT_ALWAYS(deg > 0);
   GEN_QUAD_ASSERT_ALWAYS(dim1 > 1);
   GEN_QUAD_ASSERT_ALWAYS(dim2 > 1);

   auto qs1_h1 = ComputeSimplex(deg, dim1, search_width);
   const auto &qs1 = qs1_h1.first;
   const auto &hist1 = qs1_h1.second;

   StdVector<History> hist;
   QuadSimplex qs2(deg, dim2);
   if(dim1 != dim2) {
      auto q_h2 = ComputeSimplex(deg, dim2, search_width);
      const auto &q = q_h2.first;
      const auto &his2 = q_h2.second;

      hist = dim1 > dim2 ? hist1 : his2;
      qs2.reinit_copy(q);
   } else {
      qs2.reinit_copy(qs1);
      hist = hist1;
   }

   QuadSimplexSimplex qss_init(deg, dim1, dim2, mult_nodes(qs1, qs2));
   MixedTensor(qs1, qs2, qss_init);
   auto update = NodeElimination(qss_init, search_width);
   const auto &qss_update = dynamic_cast<QuadSimplexSimplex &>(*update.first);
   hist.push_back(update.second);

   return {qss_update, hist};
}


std::pair<QuadPyramid3D, StdVector<History>> ComputePyramid3D(gq_int deg, gq_int search_width) {
   GEN_QUAD_ASSERT_ALWAYS(deg > 0);

   auto qh = ComputeCube(deg, 2, search_width);
   const auto &qc_cur = qh.first;
   auto &hist = qh.second;

   // initial guess for the pyramid
   QuadInterval qgj_next = quadrature_gauss_jacobi(deg, 0, 2);
   QuadPyramid3D qp_next(deg, mult_nodes(qgj_next, qc_cur));
   AddLinePyramid(qgj_next, qc_cur, qp_next);

   auto update = NodeElimination(qp_next, search_width);
   const auto &qp_update = dynamic_cast<QuadPyramid3D &>(*update.first);
   hist.push_back(update.second);

   return {qp_update, hist};
}


std::pair<QuadOmega2D, StdVector<History>> ComputePentagon(gq_int deg, gq_int search_width) {
   gq_int dim = 2;
   return ComputeOmega2D(CreatePentagon(), deg, search_width);
}


std::pair<QuadOmega2D, StdVector<History>> ComputeHexagon(gq_int deg, gq_int search_width) {
   return ComputeOmega2D(CreateHexagon(), deg, search_width);
}


std::pair<QuadOmega2D, StdVector<History>> ComputeOmega2D(Omega2D omega, gq_int deg, gq_int search_width) {
   gq_int dim = 2;

   QuadSimplex qs = ComputeSimplex(deg, dim, search_width).first;
   QuadOmega2D quad_omega = CreateOmegaComposite(std::move(omega), qs);
   auto update = NodeElimination(quad_omega, search_width);

   const auto &q_update = dynamic_cast<QuadOmega2D &>(*update.first);
   StdVector<History> hist;
   hist.push_back(update.second);

   return {q_update, hist};
}


}  // namespace gquad

#endif

