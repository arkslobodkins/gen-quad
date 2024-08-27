// Arkadijs Slobodkins
// Johannes Tausch
// SMU Mathematics
// Copyright 2022

#if __cplusplus < 201402L
#error requires c++14 or higher
#else

#include "../include/nonlinear_solve.hpp"

#include "../include/function.hpp"
#include "../include/util.hpp"

namespace gquad {


double LsqSolveTimer::lsq_time_total = 0.;


static void LsqUpdate(const Matrix2D& J, const Array1D& f, Array1D& dz, Array1D& dg);
static auto ExcludeIndex(gq_int index, const Array1D& x);


namespace {
struct Penalty {
   double factor;
   double distance;
};
}  // namespace


static Penalty PenaltyGradientFac(gq_int skip, const QuadConvexPolytope& q, double dz_fac,
                                  const QuadArray& dz, const QuadArray& dg);
static double PenaltyGradientFac(const QuadOmega2D& q, const QuadArray& dz, const QuadArray& dg);


static QuadArray ComputePenaltyGradient(const QuadConvexPolytope& q);
static QuadArray ComputePenaltyGradient(const QuadOmega2D& q, bool do_inner);


static std::tuple<bool, gq_int, double> PolytopeLineSearch(gq_int m, double t_prev, const Array2D& alf,
                                                           const Array2D& bet);


static Penalty OmegaLineSearch(gq_int nsteps, double tmax, const QuadOmega2D& q, const Array1D& dz);


LsqOut LeastSquaresNewton(QuadDomain& q, Basis& basis) {
   GEN_QUAD_ASSERT_ALWAYS(q.num_nodes() >= 1);

   EvalJacobian eval_jacobian(basis, orthogonal);
   EvalFunction eval_function(basis, orthogonal);

   gq_int M{basis.size()};
   gq_int N{q.size()};
   Array1D F(M);
   Matrix2D J(M, N);

   QuadArray dz(q.deg(), q.dim(), q.num_nodes());
   QuadArray dg(q.deg(), q.dim(), q.num_nodes());

   gq_int its = 0;
   gq_int max_iter = 25;

   // return if input is a satisfactory quadrature
   F = eval_function(q);  // 0th iteration
   double error_norm = function_residual(q, basis, monomial);
   if((error_norm <= constants::quad_tol) && (in_constraint(q))) {
      return LsqOut{SolFlag::found, its, error_norm};
   }

   eval_jacobian(q, J);  // 0th iteration, function for 0th iteration already computed above
   if(auto q_poly = dynamic_cast<QuadConvexPolytope*>(&q)) {
      dg = ComputePenaltyGradient(*q_poly);  // 0th iteration
   } else if(auto q_omega = dynamic_cast<QuadOmega2D*>(&q)) {
      dg = ComputePenaltyGradient(*q_omega, false);
   } else {
      GQ_THROW_RUNTIME_ERROR_MSG("Attempted domain is currently not supported by LeastSquaresNewton");
   }

   while((its < max_iter) && (error_norm > constants::quad_tol)) {
      ++its;

      util::timer t;
      LsqUpdate(J, F, dz.array(), dg.array());
      LsqSolveTimer::lsq_time_total += t.wall_time();

      double dg_fac = 0.;
      double dz_fac = 1.;
      double dist_domain = 1.;
      // for polytopes with inequality constraints
      if(auto q_poly = dynamic_cast<QuadConvexPolytope*>(&q)) {
         dz_fac = 2.;
         do {
            dz_fac *= 0.5;
            Penalty penalty = PenaltyGradientFac(-1, *q_poly, dz_fac, dz, dg);
            dg_fac = penalty.factor;
            dist_domain = penalty.distance;
         } while(dist_domain >= 0. && dz_fac > 0.05);
      }
      // for arbitrary 2D polytopes
      else if(const auto* q_omega = dynamic_cast<const QuadOmega2D*>(&q)) {
         dg_fac = PenaltyGradientFac(*q_omega, dz, dg);
         dz_fac = 2.;
         auto q_next = *q_omega;
         do {
            dz_fac *= 0.5;
            q_next.array() = q_omega->array() - dz_fac * dz.array() - dg_fac * dg.array();
            dist_domain = dist_from_constr_min(q_next);
         } while(dist_domain >= 0. && dz_fac > 0.05);
      }

      if(dist_domain > 0) {
         return {SolFlag::not_found, its, error_norm};
      }

      q.array() -= dz_fac * dz.array() + dg_fac * dg.array();

      eval_jacobian(q, J);                                       // update J
      F = eval_function(q);                                      // update F
      if(auto q_poly = dynamic_cast<QuadConvexPolytope*>(&q)) {  // update G
         dg = ComputePenaltyGradient(*q_poly);
      } else if(auto q_omega = dynamic_cast<QuadOmega2D*>(&q)) {
         dg = ComputePenaltyGradient(*q_omega, false);
      }

      error_norm = function_residual(q, basis, monomial);  // numerically more accurate
   }

   if((error_norm <= constants::quad_tol) && (in_constraint(q))) {
      return {SolFlag::found, its, error_norm};
   } else {
      return {SolFlag::not_found, its, error_norm};
   }
}


static void LsqUpdate(const Matrix2D& J, const Array1D& f, Array1D& dz, Array1D& dg) {
   GEN_QUAD_ASSERT_DEBUG(J.rows() <= J.cols());
   GEN_QUAD_ASSERT_DEBUG(dz.size() == J.cols());
   GEN_QUAD_ASSERT_DEBUG(dg.size() == dz.size());

   Eigen::HouseholderQR<Matrix2D> qr(J.transpose());
   Vector1D dzv = dz;
   Vector1D gv = dg;
   Vector1D fv = f;

   dzv = qr.transpose().solve(fv);
   dz = dzv.array();

   gv = qr.householderQ().transpose() * gv;
   gv.array()(seq(J.rows(), last)) = 0.;
   gv = qr.householderQ() * gv;

   dg -= gv.array();
}


static auto ExcludeIndex(gq_int index, const Array1D& x) {
   GEN_QUAD_ASSERT_DEBUG(index > -1 && index < x.size());
   StdVector<gq_int> indexes;
   indexes.reserve(x.size() - 1);

   for(gq_int i = 0; i < x.size(); ++i) {
      if(i != index) {
         indexes.push_back(i);
      }
   }
   return x(indexes);
}


// This routine tries to find g_fac for given dz_fac such that all nodes
// are inside the domain for qnext = q - dz_fac * dz - g_fac * dg.
// The returned value contains g_fac and the distance of the nearest node from
// the nearest boundary(including the weights). If distance < 0, all nodes are
// inside the domain, if distance > 0 then at least one node is outside of the domain.
static Penalty PenaltyGradientFac(gq_int skip, const QuadConvexPolytope& q, double dz_fac,
                                  const QuadArray& dz, const QuadArray& dg) {
   GEN_QUAD_ASSERT_DEBUG(skip >= -1 && skip <= q.num_nodes() - 1);

   const ConvexPolytope& polytope = q.get_convex_polytope();
   const Matrix2D& A = polytope.get_constraint_matrix();
   const Vector1D& b = polytope.get_constraint_vector();

   const gq_int num_constr = A.rows();
   const gq_int num_constr1 = num_constr + 1;
   const gq_int num_nodes_sk = (skip == -1) ? q.num_nodes() : q.num_nodes() - 1;

   Array1D new_node(q.dim());
   Array2D alf(num_nodes_sk, num_constr1);
   Array2D bet(num_nodes_sk, num_constr1);

   for(gq_int k = 0, kk = 0; k < q.num_nodes(); ++k) {
      if(k != skip) {
         new_node = q.node(k) - dz_fac * dz.node(k);
         alf(kk, 0) = -dg[k];
         bet(kk, 0) = -q.w(k) + dz_fac * dz[k];
         for(gq_int l = 0; l < num_constr; ++l) {
            alf(kk, l + 1) = dot_prod(A.row(l), dg.node(k));
            bet(kk, l + 1) = dot_prod(A.row(l), new_node) - b[l];
         }
         ++kk;
      }
   }

   double dz_norm = 0., g_norm = 0.;
   if(skip == -1) {
      dz_norm = norm2(dz.array());
      g_norm = norm2(dg.array());
   } else {
      dz_norm = norm2(ExcludeIndex(skip, dz.array()));
      g_norm = norm2(ExcludeIndex(skip, dg.array()));
   }
   double tmax = dz_norm / g_norm;

   bool flag;
   double t = 0.;
   gq_int m = max_index(bet);
   while(alf(m) > 0. && t < tmax) {
      std::tie(flag, m, t) = PolytopeLineSearch(m, t, alf, bet);
      if(!flag) {
         break;
      }
   }

   double distance = bet(m) - t * alf(m);
   return Penalty{t, distance};
}


static double PenaltyGradientFac(const QuadOmega2D& q, const QuadArray& dz, const QuadArray& dg) {
   auto g = ComputePenaltyGradient(q, false);
   double g_norm = norm2(g.array());
   double dz_norm = norm2(dz.array());
   double dg_norm = norm2(dg.array());

   double factor = 0.;
   if(dg_norm > 1e-5 * g_norm) {
      double t = dz_norm;
      auto qz = q;
      qz.array() -= dz.array();

      auto penalty = OmegaLineSearch(20, t, qz, dg.array() / dg_norm);
      factor = penalty.factor / dg_norm;
   }
   return factor;
}


static QuadArray ComputePenaltyGradient(const QuadConvexPolytope& q) {
   QuadArray g(q.deg(), q.dim(), q.num_nodes());

   const ConvexPolytope& polytope = q.get_convex_polytope();
   const Matrix2D& A = polytope.get_constraint_matrix();
   const Vector1D& b = polytope.get_constraint_vector();
   Vector1D denom(A.rows());

   g.weights() = -1. * norm_inf(q.weights()) / q.weights();  // contribution from weights
   for(gq_int i = 0; i < q.num_nodes(); ++i) {
      denom = b - A * q.node_v(i);
      for(gq_int j = 0; j < A.rows(); ++j) {
         for(gq_int d = 0; d < q.dim(); ++d) {
            g.node(i)(d) += A(j, d) / denom(j);
         }
      }
   }
   return g;
}


static QuadArray ComputePenaltyGradient(const QuadOmega2D& q, bool do_inner) {
   QuadArray grad(q.deg(), q.dim(), q.num_nodes());
   const Omega2D& omega = q.get_domain_self();

   for(gq_int n = 0; n < q.num_nodes(); ++n) {
      const auto xn = q.node(n);
      auto gn = grad.node(n);

      // boundary part
      for(gq_int k = 0; k < omega.edg.size(); ++k) {
         const auto& edg = omega.edg[k];
         const auto params = edg.getparams();

         if(omega.edg_on_boundary[k]) {
            const auto& v0 = edg.vert[0];
            const auto& v1 = edg.vert[1];

            const StaticArray1D<2> r{xn[0] - v0.x[0], xn[1] - v0.x[1]};
            const double r2n = std::sqrt(r[0] * r[0] + r[1] * r[1]);
            const double rt = r[0] * params.tangent[0] + r[1] * params.tangent[1];

            const StaticArray1D<2> s{xn[0] - v1.x[0], xn[1] - v1.x[1]};
            const double s2n = std::sqrt(s[0] * s[0] + s[1] * s[1]);
            const double st = s[0] * params.tangent[0] + s[1] * params.tangent[1];

            const StaticArray1D<2> res{r[0] * s2n + s[0] * r2n, r[1] * s2n + s[1] * r2n};
            const double res_dp = res[0] * res[0] + res[1] * res[1];
            if(res_dp > 1.e-20 && r2n > 1.e-12 && s2n > 1.e-12) {
               const double top = s2n - st;
               const double bot = r2n - rt;
               gn[0] += (s[0] / s2n - params.tangent[0]) / top;
               gn[1] += (s[1] / s2n - params.tangent[1]) / top;
               gn[0] -= (r[0] / r2n - params.tangent[0]) / bot;
               gn[1] -= (r[1] / r2n - params.tangent[1]) / bot;
            } else {
               gn[0] = std::numeric_limits<double>::infinity();
               gn[1] = std::numeric_limits<double>::infinity();
            }
         }
      }

      // interior part
      if(do_inner) {
         for(gq_int m = 0; m < q.num_nodes(); ++m) {
            if(n != m) {
               const auto xm = q.node(m);
               const StaticArray1D<2> r{xm[0] - xn[0], xm[1] - xn[1]};
               const double rdp = r[0] * r[0] + r[1] * r[1];
               gn[0] += r[0] / rdp;
               gn[1] += r[1] / rdp;
            }
         }
      }
   }

   return grad;
}


static std::tuple<bool, gq_int, double> PolytopeLineSearch(gq_int m, double t_prev, const Array2D& alf,
                                                           const Array2D& bet) {
   bool flag = false;
   gq_int m_new = 0;
   double t_new = std::numeric_limits<double>::infinity();

   for(gq_int k = 0; k < alf.size(); ++k) {
      if(k != m) {
         double t = (bet(k) - bet(m)) / (alf(k) - alf(m));
         if(t > t_prev && t < t_new) {
            flag = true;
            m_new = k;
            t_new = t;
         }
      }
   }

   return {flag, m_new, t_new};
}


static Penalty OmegaLineSearch(gq_int nsteps, double tmax, const QuadOmega2D& q, const Array1D& dz) {
   double t_opt = 0.0;
   double dist_opt = dist_from_boundary_min(q);
   QuadOmega2D qnext = q;
   for(gq_int k = 1; k <= nsteps; ++k) {
      double t = k * (tmax / nsteps);
      qnext.array() = q.array() - t * dz;
      double dist = dist_from_boundary_min(qnext);

      if(dist > dist_opt) {
         t_opt = t;
         dist_opt = dist;
      }
   }

   return {t_opt, dist_opt};
}


}  // namespace gquad

#endif

