// Arkadijs Slobodkins
// Johannes Tausch
// SMU Mathematics
// Copyright 2022

#if __cplusplus < 201402L
#error requires c++14 or higher
#else

#include "../include/node_elimination.hpp"

#include "../include/function.hpp"
#include "../include/nonlinear_solve.hpp"
#include "../include/util.hpp"

namespace gquad {


double PredictorTimer::predictor_time_total = 0.;


struct Distance {
   gq_int index;
   double val;
};


static void InitialSearch(QuadDomain& q, Basis& basis);


static std::pair<SolFlag, ElimData> WideLsqSearch(QuadPolytope& q, Basis& basis, gq_int search_width);
static std::pair<Matrix2D, StdVector<Distance>> Predictor(const QuadDomain& q, Basis& basis);
static Matrix2D MK_Update(const QuadDomain& q, Matrix2D J);
static void ExtractFromPredictor(gq_int array_index, const Matrix2D& Z, QuadDomain& q);


static History InitHist(const QuadDomain& q, Basis& basis);
static void IncrementHist(const ElimData& elim_data, History& h);
static void FinalizeHist(const QuadDomain& q, Basis& basis, History& h);
static void PrintElimInfo(const QuadDomain& q);
static void PrintNotFoundInfo();


static std::ostream& operator<<(std::ostream& os, const Distance& d);
static std::ostream& operator<<(std::ostream& os, const StdVector<Distance>& d);


std::pair<std::unique_ptr<QuadDomain>, History> NodeElimination(const QuadDomain& quad_init,
                                                                gq_int search_width) {
   if(!in_constraint(quad_init)) {
      GQ_THROW_RUNTIME_ERROR_MSG("Initial guess must satisfy constraints for all nodes and weights");
   }

   auto quad_new = quad_init.clone_quad_domain();
   auto basis = quad_new->create_basis();

   util::print("\neliminating nodes for " + quad_new->get_domain().domain_name());
   InitialSearch(*quad_new, *basis);
   PrintElimInfo(*quad_new);
   History h = InitHist(*quad_new, *basis);

   // run Node Elimination Algorithm. Theoretical optimum is reached when
   // (dim+1)*n_cur = num_funcs, at which point elimination is pursued no further.
   const gq_int n_opt = quad_new->num_nodes_optimal();
   gq_int n_cur = quad_new->num_nodes();
   while((n_cur > n_opt) && (n_cur >= 2)) {
      SolFlag sol_flag;
      ElimData elim_data;
      if(auto qp = dynamic_cast<QuadPolytope*>(&(*quad_new))) {
         std::tie(sol_flag, elim_data) = WideLsqSearch(*qp, *basis, search_width);
      } else {
         GQ_THROW_RUNTIME_ERROR_MSG("Attempted domain is currently not supported by NodeElimination");
      }

      if(sol_flag == SolFlag::found) {
         PrintElimInfo(*quad_new);
         IncrementHist(elim_data, h);
      } else if(sol_flag == SolFlag::not_found) {
         PrintNotFoundInfo();
         break;  // break while loop
      }

      n_cur = quad_new->num_nodes();
   }
   FinalizeHist(*quad_new, *basis, h);

   if(!in_constraint(*quad_new)) {
      GQ_THROW_RUNTIME_ERROR_MSG("UNEXPECTED RESULT: quadrature does not satisfy constraints");
   }

   util::print(function_residual(*quad_new, *basis, orthogonal),
               "final orthogonal basis residual in node_elimination");
   util::print(function_residual(*quad_new, *basis, monomial),
               "final monomial basis residual in node_elimination");

   return {std::move(quad_new), h};
}


// test accuracy of the initial quadrature
// if residual is too large, attempt to find quadrature using Newton's method.
static void InitialSearch(QuadDomain& q, Basis& basis) {
   const double tol = constants::quad_tol;  // 10^(-14);
   const double res = function_residual(q, basis, monomial);
   util::print(res, "initial monomial basis residual in node_elimination");
   if(std::fabs(res) > tol) {
      LsqOut lsq_out = LeastSquaresNewton(q, basis);
      if(lsq_out.sol_flag == SolFlag::not_found) {
         GQ_THROW_RUNTIME_ERROR_MSG("initial quadrature did not converge, bad initial guess");
      }
   }
}


static std::pair<SolFlag, ElimData> WideLsqSearch(QuadPolytope& q, Basis& basis, gq_int search_width) {
   const gq_int n_cur{q.num_nodes()};
   const gq_int max_fail_wide{20};
   const gq_int max_fail = (q.dim() == 2) ? n_cur : 100;
   const gq_int max_success = (q.efficiency() > 0.5) ? std::min(search_width, n_cur) : 1;

   Matrix2D Z{};
   StdVector<Distance> dst{};
   std::tie(Z, dst) = Predictor(q, basis);

   StdVector<gq_int> good_indexes{};
   StdVector<LsqOut> good_lsq{};
   StdVector<std::unique_ptr<QuadPolytope>> good_quad{};

   auto q_temp = q.clone_quad_polytope();
   q_temp->resize(n_cur - 1);

   gq_int good_count = 0;
   gq_int fail_count = 0;
   for(gq_int i = 0; i < n_cur; ++i) {
      if(dst[i].val > constants::quad_huge) {
         continue;  // catches large dz and not in domain
      }
      ExtractFromPredictor(dst[i].index, Z, *q_temp);
      auto lsq_cur = LeastSquaresNewton(*q_temp, basis);

      if(lsq_cur.sol_flag == SolFlag::found) {
         GEN_QUAD_ASSERT_DEBUG(in_constraint(*q_temp));

         // discard really small weights to avoid complications; neither good_count nor fail_count are
         // incremented
         if(q_temp->weights().minCoeff() < 1.e-12) {
            util::print("Warning: Found solution containing one or more small weights");
            continue;
         }

         good_indexes.push_back(i);
         good_lsq.push_back(lsq_cur);
         good_quad.push_back(q_temp->clone_quad_polytope());
         if(++good_count >= max_success) {
            break;
         }
      } else {
         // break due to max_fail_wide only if at least one quadrature has succeeded,
         // otherwise continue looping over more iterations specified by max_fail
         if((++fail_count >= max_fail_wide && good_count > 0) || (fail_count >= max_fail)) {
            break;
         }
      }
   }

   if(good_count == 0) {
      return {SolFlag::not_found, ElimData{}};
   }

   Array1D dist(good_count);
   for(gq_int i = 0; i < good_count; ++i) {
      dist[i] = -dist_from_constr_min(*good_quad[i]);
      GEN_QUAD_ASSERT_DEBUG(dist[i] >= 0.);
   }
   gq_int max_dist_index = max_index(dist);
   q.resize_and_assign(*good_quad[max_dist_index]);

   ElimData data;
   data.nodes_tot = q.num_nodes();
   data.success_node = good_indexes.at(max_dist_index);
   data.success_its = good_lsq[max_dist_index].its;
   data.num_solutions = good_count;
   data.num_fails = fail_count;

   return {SolFlag::found, data};
}


// Calculates all predictors obtained by setting each weight to zero.
// Returns a matrix that contains in its rows the z-vectors of predicted values
// and Distance vector which sorts the predictors according to their distance.
// Here, nodes which are not in the domain have a large distance, appear at the end
// and their z-vectors are not defined.
static std::pair<Matrix2D, StdVector<Distance>> Predictor(const QuadDomain& q, Basis& basis) {
   const gq_int M = basis.size();
   const gq_int N = q.size();

   Matrix2D J(M, N);
   EvalJacobian eval_jacobian(basis, orthogonal);
   eval_jacobian(q, J);

   util::timer t;  // time from here since Jacobian is timed separately

   Matrix2D Q2 = MK_Update(q, std::move(J));
   Matrix2D MK = Q2(all, seqN(0, q.num_nodes()));
   Matrix2D MP = Q2.transpose() * MK;

   Matrix2D Z(q.size(), q.num_nodes());
   StdVector<Distance> dst(q.num_nodes());

   auto qz = q.clone_quad_domain();

   for(gq_int k = 0; k < q.num_nodes(); ++k) {
      const double cn2 = norm2(MK.col(k));
      const double z_fac = q.w(k) / (cn2 * cn2);
      const auto dz = MP.col(k) * z_fac;

      qz->array() = q.array() - dz.array();
      if(in_constraint(*qz)) {
         Z.col(k) = qz->array();
         dst[k].index = k;
         dst[k].val = norm2(dz.array());
      } else {
         dst[k].val = std::numeric_limits<double>::infinity();
      }
   }

   std::stable_sort(dst.begin(), dst.end(), [](Distance a, Distance b) { return a.val < b.val; });

   PredictorTimer::predictor_time_total += t.wall_time();
   return {Z, dst};
}


static Matrix2D MK_Update(const QuadDomain& q, Matrix2D J) {
   GEN_QUAD_ASSERT_ALWAYS(J.cols() > J.rows());

   gq_int M{J.rows()};
   gq_int N{J.cols()};

   Eigen::HouseholderQR<Matrix2D> qr(J.transpose());  // "N x N"
   Matrix2D Q2(N, N - M);
   Q2(seq(M, last), all) = Eigen::MatrixXd::Identity(N - M, N - M);
   Q2 = (qr.householderQ() * Q2).transpose();  // (N - M) x N

   return Q2;
}


static void ExtractFromPredictor(const gq_int array_index, const Matrix2D& Z, QuadDomain& q) {
   GEN_QUAD_ASSERT_DEBUG(array_index > -1 && array_index < Z.cols());
   GEN_QUAD_ASSERT_DEBUG((q.num_nodes() + 1) * (q.dim() + 1)
                         == Z.rows());  // q must contain one fewer nodes than previous quadrature

   gq_int dim{q.dim()};
   gq_int num_nodes{q.num_nodes()};

   for(gq_int j = 0; j < array_index; ++j) {
      q.w(j) = Z(j, array_index);
   }

   for(gq_int j = array_index + 1; j < num_nodes + 1; ++j) {
      q.w(j - 1) = Z(j, array_index);
   }

   for(gq_int j = 0; j < array_index; ++j) {
      for(gq_int d = 0; d < dim; ++d) {
         q.node(j)(d) = Z((num_nodes + 1) + j * dim + d, array_index);
      }
   }

   for(gq_int j = array_index + 1; j < num_nodes + 1; ++j) {
      for(gq_int d = 0; d < dim; ++d) {
         q.node(j - 1)(d) = Z((num_nodes + 1) + j * dim + d, array_index);
      }
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
static History InitHist(const QuadDomain& q, Basis& basis) {
   History h{};
   h.deg = q.deg();
   h.dim = q.dim();
   h.nodes_initial = q.num_nodes();
   h.nodes_optimal = q.num_nodes_optimal();
   h.num_funcs = basis.size();
   h.total_elims = 0;
   h.domain_name = q.get_domain().domain_name();

   return h;
}


static void IncrementHist(const ElimData& elim_data, History& h) {
   h.hlist.push_back(elim_data);
   ++h.total_elims;
}


static void FinalizeHist(const QuadDomain& q, Basis& basis, History& h) {
   h.nodes_final = q.num_nodes();
   h.efficiency = double(h.nodes_optimal) / h.nodes_final;
   h.res = function_residual(q, basis, monomial);
}


static void PrintElimInfo(const QuadDomain& q) {
   gq_int dim{q.dim()};
   gq_int num_nodes{q.num_nodes()};
   gq_int n_opt{q.num_nodes_optimal()};
   double efficiency{q.efficiency()};
   std::printf("dimension = %ld, current number of nodes = %ld, optimal = %ld, efficiency = %f\n",
               dim,
               num_nodes,
               n_opt,
               efficiency);
}


static void PrintNotFoundInfo() {
   util::print("Last iteration did not converge");
}


////////////////////////////////////////////////////////////////////////////////////////////////////
static std::ostream& operator<<(std::ostream& os, const Distance& d) {
   std::cout << "index: " << d.index << " value: " << d.val << std::endl;
   return os;
}


static std::ostream& operator<<(std::ostream& os, const StdVector<Distance>& dist) {
   os << std::scientific << std::showpoint << std::setprecision(3);

   using namespace util;
   for(gq_int i = 0; i < dist.size(); ++i) {
      os << index_of("dist", i) + ":" << smart_spaces(gq_int(dist.size() - 1), i)
         << "index: " << dist[i].index << smart_spaces(gq_int(dist.size() - 1), dist[i].index)
         << " value: " << dist[i].val << std::endl;
   }
   return os;
}


}  // namespace gquad

#endif

