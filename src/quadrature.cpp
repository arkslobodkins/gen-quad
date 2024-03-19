// Arkadijs Slobodkins
// Johannes Tausch
// SMU Mathematics
// Copyright 2022

#if __cplusplus < 201402L
#error requires c++14 or higher

#else

#include "../include/quadrature.hpp"

#include "../alglib/src/integration.h"
#include "../include/math_func.hpp"
#include "../include/quad_tensor.hpp"
#include "../include/util.hpp"

namespace gquad {

static void assert_sizes(gq_int deg, gq_int dim, gq_int num_nodes);
static void assert_sizes(const QuadArray& q1, const QuadArray& q2);

static double rel_error(double approx, long double exact);

static long double expIntegralNDimCube(gq_int dim);
static long double expIntegralNDimSimplex(gq_int dim);

QuadArray::QuadArray(gq_int deg, gq_int dim, gq_int num_nodes)
    : m_array(num_nodes * (dim + 1)),
      m_deg{deg},
      m_dim{dim},
      m_num_nodes{num_nodes} {
   assert_sizes(m_deg, m_dim, m_num_nodes);
}

QuadArray::QuadArray(const QuadArray& q) : QuadArray{q.deg(), q.dim(), q.num_nodes()} {
   m_array = q.m_array;
}

QuadArray::QuadArray(QuadArray&& q)
    : m_array{std::move(q.m_array)},
      m_deg{q.m_deg},
      m_dim{q.m_dim},
      m_num_nodes{q.m_num_nodes} {
   q.m_array = {};
   q.m_deg = 0;
   q.m_dim = 0;
   q.m_num_nodes = 0;
}

QuadArray& QuadArray::operator=(const QuadArray& q) & {
   if(this != &q) {
      assert_sizes(*this, q);
      GEN_QUAD_ASSERT_DEBUG(m_array.size() == q.m_array.size());
      m_array = q.m_array;
   }
   return *this;
}

QuadArray& QuadArray::operator=(QuadArray&& q) & {
   if(this != &q) {
      assert_sizes(*this, q);
      GEN_QUAD_ASSERT_DEBUG(m_array.size() == q.m_array.size());

      m_array = std::move(q.m_array);
      q.m_array = {};
      q.m_deg = 0;
      q.m_dim = 0;
      q.m_num_nodes = 0;
   }

   return *this;
}

void QuadArray::resize(gq_int new_num_nodes) {
   GEN_QUAD_ASSERT_DEBUG(new_num_nodes > -1);
   m_array.resize(new_num_nodes * (m_dim + 1));
   m_num_nodes = new_num_nodes;
}

void QuadArray::resize_and_assign(const QuadArray& q) {
   this->resize(q.num_nodes());
   *this = q;
}

void QuadArray::reinit(gq_int new_dim, gq_int new_num_nodes) {
   GEN_QUAD_ASSERT_DEBUG(new_dim > 0);
   GEN_QUAD_ASSERT_DEBUG(new_num_nodes > -1);
   m_array.resize(new_num_nodes * (new_dim + 1));
   m_dim = new_dim;
   m_num_nodes = new_num_nodes;
}

std::ostream& operator<<(std::ostream& os, const QuadArray& q) {
   os << "deg: " << q.deg() << std::endl;
   os << "dim: " << q.dim() << std::endl;
   os << "nodes: " << q.num_nodes() << std::endl;
   os << std::scientific << std::showpoint << std::setprecision(16);

   using namespace util;
   for(gq_int i = 0; i < q.num_nodes(); ++i) {
      os << index_of("q", i) + ":" << smart_spaces(q.num_nodes() - 1, i) << q.w(i) << n_spaces(2);
      for(gq_int d = 0; d < q.dim(); ++d) {
         os << q.node(i)(d) << n_spaces(2);
      }
      os << std::endl;
   }
   return os;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// abtract class
QuadDomain::QuadDomain(gq_int deg, gq_int dim, gq_int num_nodes) : QuadArray(deg, dim, num_nodes) {
}

std::unique_ptr<Basis> QuadDomain::create_basis() const {
   GQ_THROW_RUNTIME_ERROR();
   return nullptr;
}

double QuadDomain::efficiency() const {
   return double(num_nodes_optimal()) / num_nodes();
}

QuadDomain& QuadDomain::resize_and_assign(const QuadDomain& q) {
   QuadArray::resize_and_assign(q);
   return *this;
}

bool in_domain_elem(const QuadDomain& q, gq_int i) {
   const auto& domain = q.get_domain();
   return domain.in_domain(q.node(i));
}

bool in_domain(const QuadDomain& q) {
   GEN_QUAD_ASSERT_DEBUG(q.num_nodes() > 0);
   for(gq_int i = 0; i < q.num_nodes(); ++i) {
      if(!in_domain_elem(q, i)) {
         return false;
      }
   }
   return true;
}

bool pos_weights(const QuadDomain& q) {
   GEN_QUAD_ASSERT_DEBUG(q.num_nodes() > 0);
   for(gq_int i = 0; i < q.num_nodes(); ++i) {
      if(q.w(i) < 0.) {
         return false;
      }
   }
   return true;
}

bool in_constraint(const QuadDomain& q) {
   return in_domain(q) && pos_weights(q);
}

std::ostream& operator<<(std::ostream& os, const QuadDomain& q) {
   os << QuadArray(q);
   return os;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// abstract class
QuadPolytope::QuadPolytope(gq_int deg, gq_int dim, gq_int num_nodes) : QuadDomain(deg, dim, num_nodes) {
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// abstract class
QuadIdealPolytope::QuadIdealPolytope(gq_int deg, gq_int dim, gq_int num_nodes)
    : QuadPolytope(deg, dim, num_nodes) {
}

double dist_from_boundary_min(const QuadPolytope& q) {
   if(!in_domain(q)) {
      return std::numeric_limits<double>::infinity();
   }

   Vector1D distances(q.num_nodes());
   const auto& polytope = q.get_polytope();

   for(gq_int i = 0; i < q.num_nodes(); ++i) {
      distances[i] = polytope.dist_from_boundary(q.node(i));
   }
   return distances.maxCoeff();
}

double dist_from_constr_min(const QuadPolytope& q) {
   if(!in_constraint(q)) {
      return std::numeric_limits<double>::infinity();
   }

   double weight_dist
       = -q.weights().minCoeff() * q.num_nodes() * q.dim() * q.dim();  // scale the smallest weight
   double domain_dist = dist_from_boundary_min(q);
   return std::max(weight_dist,
                   domain_dist);  // max of two negative values, i.e. min of their absolute values
}

////////////////////////////////////////////////////////////////////////////////////////////////////
QuadInterval::QuadInterval(gq_int deg, gq_int num_nodes) : QuadIdealPolytope(deg, 1, num_nodes) {
}

std::unique_ptr<QuadDomain> QuadInterval::clone_quad_domain() const {
   return std::make_unique<QuadInterval>(QuadInterval(*this));
}

std::unique_ptr<QuadPolytope> QuadInterval::clone_quad_polytope() const {
   return std::make_unique<QuadInterval>(QuadInterval(*this));
}

const Domain& QuadInterval::get_domain() const {
   return interval;
}

const Polytope& QuadInterval::get_polytope() const {
   return interval;
}

const IdealPolytope& QuadInterval::get_ideal_polytope() const {
   return interval;
}

double QuadInterval::relative_exponential_residual() const {
   return rel_error(quad_exp_approx(*this), expIntegralNDimCube(1));
}

std::string QuadInterval::quad_file_name() const {
   using std::string;
   using std::to_string;
   return string("quad_") + interval.domain_name() + string("_deg") + to_string(this->deg()) + ".txt";
}

std::string QuadInterval::quad_dims_file_name() const {
   return interval.domain_name() + ".txt";
}

gq_int QuadInterval::num_nodes_optimal() const {
   return 1 + deg() / 2;
}

QuadInterval quadrature_gauss_legendre(gq_int deg) {
   return quadrature_gauss_jacobi(deg, 0, 0);
}

QuadInterval quadrature_gauss_jacobi(gq_int deg, gq_int alpha, gq_int beta) {
   gq_int n{1 + deg / 2};
   QuadInterval q(deg, n);

   alglib::ae_int_t info;
   alglib::real_1d_array x, w;
   x.setlength(n);
   w.setlength(n);

   alglib::gqgenerategaussjacobi(alglib::ae_int_t(n), double(alpha), double(beta), info, x, w);
   assert(info == 1);

   double c = std::pow(2., double(alpha + beta + 1));
   for(alglib::ae_int_t i = 0; i < n; ++i) {
      q.w(i) = w[i] / c;
      q.node(i)(0) = (x[i] + 1.) / 2.;
   }

   return q;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
QuadCube::QuadCube(gq_int deg, gq_int dim, gq_int num_nodes)
    : QuadIdealPolytope(deg, dim, num_nodes),
      cube(dim) {
   GEN_QUAD_ASSERT_DEBUG(dim >= 2);
}

std::unique_ptr<QuadDomain> QuadCube::clone_quad_domain() const {
   return std::make_unique<QuadCube>(QuadCube(*this));
}

std::unique_ptr<QuadPolytope> QuadCube::clone_quad_polytope() const {
   return std::make_unique<QuadCube>(QuadCube(*this));
}

const Domain& QuadCube::get_domain() const {
   return cube;
}

const Polytope& QuadCube::get_polytope() const {
   return cube;
}

const IdealPolytope& QuadCube::get_ideal_polytope() const {
   return cube;
}

std::unique_ptr<Basis> QuadCube::create_basis() const {
   return std::make_unique<CubeBasis>(CubeBasis{this->deg(), this->dim()});
}

double QuadCube::relative_exponential_residual() const {
   return rel_error(quad_exp_approx(*this), expIntegralNDimCube(this->dim()));
}

std::string QuadCube::quad_file_name() const {
   using std::string;
   using std::to_string;
   return string("quad_") + cube.domain_name() + string("_deg") + to_string(this->deg()) + string("_dim")
        + to_string(this->dim()) + ".txt";
}

std::string QuadCube::quad_dims_file_name() const {
   using std::string;
   using std::to_string;
   return cube.domain_name() + string("_dim") + to_string(this->dim()) + ".txt";
}

QuadCube& QuadCube::reinit(gq_int dim, gq_int num_nodes) {
   GEN_QUAD_ASSERT_DEBUG(dim > 1);
   QuadArray::reinit(dim, num_nodes);
   cube.reinit(dim);
   return *this;
}

QuadCube& QuadCube::reinit_copy(const QuadCube& qc) {
   reinit(qc.dim(), qc.num_nodes());
   *this = qc;
   return *this;
}

gq_int QuadCube::num_nodes_optimal() const {
   gq_int basis_size = CubeBasisSize(deg(), dim());
   return (basis_size + dim()) / (dim() + 1);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
QuadSimplex::QuadSimplex(gq_int deg, gq_int dim, gq_int num_nodes)
    : QuadIdealPolytope(deg, dim, num_nodes),
      simplex(dim) {
   GEN_QUAD_ASSERT_DEBUG(dim >= 2);
}

std::unique_ptr<QuadDomain> QuadSimplex::clone_quad_domain() const {
   return std::make_unique<QuadSimplex>(QuadSimplex(*this));
}

std::unique_ptr<QuadPolytope> QuadSimplex::clone_quad_polytope() const {
   return std::make_unique<QuadSimplex>(QuadSimplex(*this));
}

const Domain& QuadSimplex::get_domain() const {
   return simplex;
}

const Polytope& QuadSimplex::get_polytope() const {
   return simplex;
}

const IdealPolytope& QuadSimplex::get_ideal_polytope() const {
   return simplex;
}

std::unique_ptr<Basis> QuadSimplex::create_basis() const {
   return std::make_unique<SimplexBasis>(SimplexBasis{this->deg(), this->dim()});
}

double QuadSimplex::relative_exponential_residual() const {
   return rel_error(quad_exp_approx(*this), expIntegralNDimSimplex(this->dim()));
}

std::string QuadSimplex::quad_file_name() const {
   using std::string;
   using std::to_string;
   return string("quad_") + simplex.domain_name() + string("_deg") + to_string(this->deg()) + string("_dim")
        + to_string(this->dim()) + ".txt";
}

std::string QuadSimplex::quad_dims_file_name() const {
   using std::string;
   using std::to_string;
   return simplex.domain_name() + string("_dim") + to_string(this->dim()) + ".txt";
}

QuadSimplex& QuadSimplex::reinit(gq_int dim, gq_int num_nodes) {
   GEN_QUAD_ASSERT_DEBUG(dim > 1);
   QuadArray::reinit(dim, num_nodes);
   simplex.reinit(dim);
   return *this;
}

QuadSimplex& QuadSimplex::reinit_copy(const QuadSimplex& qs) {
   reinit(qs.dim(), qs.num_nodes());
   *this = qs;
   return *this;
}

gq_int QuadSimplex::num_nodes_optimal() const {
   gq_int basis_size = SimplexBasisSize(deg(), dim());
   return (basis_size + dim()) / (dim() + 1);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
QuadCubeSimplex::QuadCubeSimplex(gq_int deg, gq_int dim1, gq_int dim2, gq_int num_nodes)
    : QuadIdealPolytope(deg, dim1 + dim2, num_nodes),
      cubesimplex({dim1, dim2}) {
   GEN_QUAD_ASSERT_DEBUG(dim1 >= 1);
   GEN_QUAD_ASSERT_DEBUG(dim2 >= 2);
}

std::unique_ptr<QuadDomain> QuadCubeSimplex::clone_quad_domain() const {
   return std::make_unique<QuadCubeSimplex>(QuadCubeSimplex(*this));
}

std::unique_ptr<QuadPolytope> QuadCubeSimplex::clone_quad_polytope() const {
   return std::make_unique<QuadCubeSimplex>(QuadCubeSimplex(*this));
}

const Domain& QuadCubeSimplex::get_domain() const {
   return cubesimplex;
}

const Polytope& QuadCubeSimplex::get_polytope() const {
   return cubesimplex;
}

const IdealPolytope& QuadCubeSimplex::get_ideal_polytope() const {
   return cubesimplex;
}

std::unique_ptr<Basis> QuadCubeSimplex::create_basis() const {
   return std::make_unique<CubeSimplexBasis>(CubeSimplexBasis{this->deg(), this->dim1(), this->dim2()});
}

double QuadCubeSimplex::relative_exponential_residual() const {
   return rel_error(quad_exp_approx(*this),
                    expIntegralNDimCube(this->dim1()) * expIntegralNDimSimplex(this->dim2()));
}

std::string QuadCubeSimplex::quad_file_name() const {
   using std::string;
   using std::to_string;
   return string("quad_") + cubesimplex.domain_name() + string("_deg") + to_string(this->deg())
        + string("_dims_") + to_string(this->dim1()) + string("_") + to_string(this->dim2()) + ".txt";
}

std::string QuadCubeSimplex::quad_dims_file_name() const {
   using std::string;
   using std::to_string;
   return cubesimplex.domain_name() + string("_dims_") + to_string(this->dim1()) + string("_")
        + to_string(this->dim2()) + ".txt";
}

QuadCubeSimplex& QuadCubeSimplex::reinit(gq_int dim1, gq_int dim2, gq_int num_nodes) {
   GEN_QUAD_ASSERT_DEBUG(dim1 > 0);
   GEN_QUAD_ASSERT_DEBUG(dim2 > 1);
   QuadArray::reinit(dim1 + dim2, num_nodes);
   cubesimplex.reinit(dim1, dim2);
   return *this;
}

QuadCubeSimplex& QuadCubeSimplex::reinit_copy(const QuadCubeSimplex& qcs) {
   reinit(qcs.dim1(), qcs.dim2(), qcs.num_nodes());
   *this = qcs;
   return *this;
}

gq_int QuadCubeSimplex::num_nodes_optimal() const {
   gq_int basis_size = CubeSimplexBasisSize(deg(), dim());
   return (basis_size + dim()) / (dim() + 1);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
QuadSimplexSimplex::QuadSimplexSimplex(gq_int deg, gq_int dim1, gq_int dim2, gq_int num_nodes)
    : QuadIdealPolytope(deg, dim1 + dim2, num_nodes),
      simplexsimplex({dim1, dim2}) {
   GEN_QUAD_ASSERT_DEBUG(dim1 >= 2);
   GEN_QUAD_ASSERT_DEBUG(dim2 >= 2);
}

std::unique_ptr<QuadDomain> QuadSimplexSimplex::clone_quad_domain() const {
   return std::make_unique<QuadSimplexSimplex>(QuadSimplexSimplex(*this));
}

std::unique_ptr<QuadPolytope> QuadSimplexSimplex::clone_quad_polytope() const {
   return std::make_unique<QuadSimplexSimplex>(QuadSimplexSimplex(*this));
}

const Domain& QuadSimplexSimplex::get_domain() const {
   return simplexsimplex;
}

const Polytope& QuadSimplexSimplex::get_polytope() const {
   return simplexsimplex;
}

const IdealPolytope& QuadSimplexSimplex::get_ideal_polytope() const {
   return simplexsimplex;
}

std::unique_ptr<Basis> QuadSimplexSimplex::create_basis() const {
   return std::make_unique<SimplexSimplexBasis>(SimplexSimplexBasis{this->deg(), this->dim1(), this->dim2()});
}

double QuadSimplexSimplex::relative_exponential_residual() const {
   return rel_error(quad_exp_approx(*this),
                    expIntegralNDimSimplex(this->dim1()) * expIntegralNDimSimplex(this->dim2()));
}

std::string QuadSimplexSimplex::quad_file_name() const {
   using std::string;
   using std::to_string;
   return string("quad_") + simplexsimplex.domain_name() + string("_deg") + to_string(this->deg())
        + string("_dims_") + to_string(this->dim1()) + string("_") + to_string(this->dim2()) + ".txt";
}

std::string QuadSimplexSimplex::quad_dims_file_name() const {
   using std::string;
   using std::to_string;
   return simplexsimplex.domain_name() + string("_dims_") + to_string(this->dim1()) + string("_")
        + to_string(this->dim2()) + ".txt";
}

QuadSimplexSimplex& QuadSimplexSimplex::reinit(gq_int dim1, gq_int dim2, gq_int num_nodes) {
   GEN_QUAD_ASSERT_DEBUG(dim1 > 1);
   GEN_QUAD_ASSERT_DEBUG(dim2 > 1);
   QuadArray::reinit(dim1 + dim2, num_nodes);
   simplexsimplex.reinit(dim1, dim2);
   return *this;
}

QuadSimplexSimplex& QuadSimplexSimplex::reinit_copy(const QuadSimplexSimplex& qss) {
   reinit(qss.dim1(), qss.dim2(), qss.num_nodes());
   *this = qss;
   return *this;
}

gq_int QuadSimplexSimplex::num_nodes_optimal() const {
   gq_int basis_size = SimplexSimplexBasisSize(deg(), dim());
   return (basis_size + dim()) / (dim() + 1);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
QuadPyramid3D::QuadPyramid3D(gq_int deg, gq_int num_nodes)
    : QuadIdealPolytope(deg, 3, num_nodes),
      pyramid3D{} {
}

std::unique_ptr<QuadDomain> QuadPyramid3D::clone_quad_domain() const {
   return std::make_unique<QuadPyramid3D>(QuadPyramid3D(*this));
}

std::unique_ptr<QuadPolytope> QuadPyramid3D::clone_quad_polytope() const {
   return std::make_unique<QuadPyramid3D>(QuadPyramid3D(*this));
}

const Domain& QuadPyramid3D::get_domain() const {
   return pyramid3D;
}

const Polytope& QuadPyramid3D::get_polytope() const {
   return pyramid3D;
}

const IdealPolytope& QuadPyramid3D::get_ideal_polytope() const {
   return pyramid3D;
}

std::unique_ptr<Basis> QuadPyramid3D::create_basis() const {
   return std::make_unique<PyramidBasis3D>(PyramidBasis3D{this->deg()});
}

double QuadPyramid3D::relative_exponential_residual() const {
   return rel_error(quad_exp_approx(*this), expIntegralNDimCube(3) / 3.L);
}

std::string QuadPyramid3D::quad_file_name() const {
   using std::string;
   using std::to_string;
   return string("quad_") + pyramid3D.domain_name() + string("_deg") + to_string(this->deg()) + ".txt";
}

std::string QuadPyramid3D::quad_dims_file_name() const {
   return pyramid3D.domain_name() + ".txt";
}

gq_int QuadPyramid3D::num_nodes_optimal() const {
   gq_int basis_size = Pyramid3DBasisSize(deg());
   return (basis_size + dim()) / (dim() + 1);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
QuadOmega2D::QuadOmega2D(Omega2D omega_, gq_int deg, gq_int num_nodes)
    : QuadPolytope(deg, 2, num_nodes),
      omega{std::move(omega_)} {
}

std::unique_ptr<QuadDomain> QuadOmega2D::clone_quad_domain() const {
   return std::make_unique<QuadOmega2D>(QuadOmega2D(*this));
}

std::unique_ptr<QuadPolytope> QuadOmega2D::clone_quad_polytope() const {
   return std::make_unique<QuadOmega2D>(QuadOmega2D(*this));
}

const Domain& QuadOmega2D::get_domain() const {
   return omega;
}

const Polytope& QuadOmega2D::get_polytope() const {
   return omega;
}

const Omega2D& QuadOmega2D::get_domain_loc() const {
   return omega;
}

std::unique_ptr<Basis> QuadOmega2D::create_basis() const {
   return std::make_unique<OmegaBasis2D>(OmegaBasis2D{omega, this->deg()});
}

std::unique_ptr<OmegaBasis2D> QuadOmega2D::create_basis_loc() const {
   return std::make_unique<OmegaBasis2D>(OmegaBasis2D{omega, this->deg()});
}

double QuadOmega2D::relative_exponential_residual() const {
   QuadOmega2D qe = CreateOmegaComposite(omega, GaussTensorSimplex(2 * deg() + 1, 2));
   return rel_error(quad_exp_approx(*this), quad_exp_approx(qe));
}

std::string QuadOmega2D::quad_file_name() const {
   using std::string;
   using std::to_string;
   return string("quad_") + omega.domain_name() + string("_deg") + to_string(this->deg()) + ".txt";
}

std::string QuadOmega2D::quad_dims_file_name() const {
   return omega.domain_name() + ".txt";
}

gq_int QuadOmega2D::num_nodes_optimal() const {
   gq_int basis_size = Omega2DBasisSize(deg());
   return (basis_size + dim()) / (dim() + 1);
}

QuadOmega2D CreateOmegaComposite(Omega2D omega, const QuadSimplex& qs) {
   GEN_QUAD_ASSERT_DEBUG(qs.dim() == 2);

   const auto nt = omega.triang.size();
   const auto triang = omega.triang;
   QuadOmega2D q(std::move(omega), qs.deg(), qs.num_nodes() * nt);

   gq_int count = 0;
   for(gq_int n = 0; n < nt; ++n) {
      const auto& tr = triang[n];
      for(gq_int k = 0; k < qs.num_nodes(); ++k) {
         map_from_unit(tr, qs.node(k), q.node(count));
         q.w(count++) = qs.w(k) * tr.jacobian();
      }
   }

   return q;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
static void assert_sizes(gq_int deg, gq_int dim, gq_int num_nodes) {
   GEN_QUAD_ASSERT_DEBUG(deg > 0);
   GEN_QUAD_ASSERT_DEBUG(dim > 0);
   GEN_QUAD_ASSERT_DEBUG(num_nodes >= 0);
}

static void assert_sizes(const QuadArray& q1, const QuadArray& q2) {
   GEN_QUAD_ASSERT_DEBUG(q1.deg() == q2.deg());
   GEN_QUAD_ASSERT_DEBUG(q1.dim() == q2.dim());
   GEN_QUAD_ASSERT_DEBUG(q1.num_nodes() == q2.num_nodes());
}

static double rel_error(double approx, long double exact) {
   return static_cast<double>(::fabsl(exact - approx) / exact);
}

static long double expIntegralNDimCube(gq_int dim) {
   long double expI = 1.L;
   for(gq_int i = 0; i < dim; ++i) {
      expI *= math::expIntegral1D(1.L);
   }
   return expI;
}

static long double expIntegralNDimSimplex(gq_int dim) {
   long double expI = 0.L;
   gq_int sign = 1;
   for(gq_int i = dim; i >= 0; --i) {
      expI += (long double)math::binomial(i, dim) * std::exp(((long double)i)) * (long double)sign;
      sign *= -1;
   }
   expI /= (long double)math::factorial(dim);

   return expI;
}

double quad_exp_approx(const QuadDomain& q) {
   Array1D eval(q.num_nodes());
   for(gq_int i = 0; i < q.num_nodes(); ++i) {
      gq_int offset = q.dim() * i;
      eval[i] = math::expNDim(q.dim(), q.x_ptr() + offset);
   }

   return stable_dot_prod(eval, q.weights());
}

}  // namespace gquad

#endif

